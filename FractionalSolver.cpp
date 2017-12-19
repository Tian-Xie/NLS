/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "FractionalSolver.h"
#include <assert.h>
#include <algorithm>

const double FractionalSolverEpsilon = 1e-8;
const double FractionalSolverOptimalityEpsilon = 1e-8;
const double FractionalSolverInfinity = 1e100;

double DPOracle(int m, int* n, int C, double** A, double** B, double a0, double b0, double z, double** dp, int** prev)
{
	for (int j = 0; j <= C; j ++)
		dp[0][j] = a0 - b0 * z;

	if (prev != NULL)
	{
		for (int i = 1; i <= m; i ++)
			for (int j = 0; j <= C; j ++)
			{
				dp[i][j] = -FractionalSolverInfinity;
				for (int d = 0; d <= j && d <= n[i - 1]; d ++)
				{
					double cur_value = dp[i - 1][j - d] + A[i - 1][d] - B[i - 1][d] * z;
					if (prev[i][j] == -1 || cur_value > dp[i][j])
					{
						dp[i][j] = cur_value;
						prev[i][j] = d;
					}
				}
			}
	}
	else // an optimization if we do not care about the trajectory
	{
		for (int i = 1; i <= m; i ++)
			for (int j = 0; j <= C; j ++)
			{
				dp[i][j] = -FractionalSolverInfinity;
				for (int d = 0; d <= j && d <= n[i - 1]; d ++)
				{
					double cur_value = dp[i - 1][j - d] + A[i - 1][d] - B[i - 1][d] * z;
					if (cur_value > dp[i][j])
						dp[i][j] = cur_value;
				}
			}
	}
	return dp[m][C];
}

#ifdef _USE_CPLEX
#include <ilcplex\ilocplex.h>
#endif

void MaximizeFraction(int m, int* n, int C, double** A, double** B, double a0, double b0, double& z, int* q)
{
	int totElement = 0;
	int n_max = 0;
	for (int i = 0; i < m; i ++)
	{
		totElement += n[i];
		n_max = std::max(n_max, n[i]);
	}
	C = std::max(C, 0);
	C = std::min(C, totElement);

	// Set a upper bound and lower bound
	double zL = std::min(a0, 0.), zR = std::max(a0, 0.);
	for (int i = 0; i < m; i ++)
		for (int j = 0; j < n[i]; j ++)
		{
			zL += std::min(A[i][j], 0.);
			zR += std::max(A[i][j], 0.);
		}
	zL /= b0;
	zR /= b0;
	zL -= 1.0;
	zR += 1.0;

	double* last_a = new double[C + 1];
	double* last_b = new double[C + 1];
	double* cur_a = new double[C + 1];
	double* cur_b = new double[C + 1];
	for (int j = 0; j <= C; j ++)
	{
		last_a[j] = a0;
		last_b[j] = b0;
	}

	double** dp = new double*[m + 1];
	for (int i = 0; i <= m; i ++)
		dp[i] = new double[C + 1];
	int** prev = new int*[m + 1];
	for (int i = 0; i <= m; i ++)
		prev[i] = new int[C + 1];

	double* break_point = new double[n_max * (C + 1) + 2];
	double* cur_alpha = new double[n_max + 1];
	double* cur_beta = new double[n_max + 1];
	int* cur_order = new int[n_max + 1];
	double* cur_range = new double[n_max + 1];

	for (int i = 0; i < m; i ++)
	{
		int nbreak_point = 0;
		for (int j = 0; j <= C; j ++)
		{
			int cur_n = 0;
			for (int d = 0; d <= j && d <= n[i]; d ++)
			{
				cur_alpha[cur_n] = last_a[j - d] + A[i][d];
				cur_beta[cur_n] = last_b[j - d] + B[i][d];
				cur_order[cur_n] = cur_n;
				cur_n ++;
			}
			std::sort(cur_order, cur_order + cur_n, [&](const int& a, const int& b) -> bool
			{
				if (fabs(cur_beta[a] - cur_beta[b]) < FractionalSolverEpsilon)
					return cur_alpha[a] > cur_alpha[b];
				return cur_beta[a] > cur_beta[b];
			});

			int cur_filled = 0;
			for (int k = 0; k < cur_n; k ++)
			{
				while (cur_filled >= 1)
				{
					if (cur_alpha[cur_order[cur_filled - 1]] - zR * cur_beta[cur_order[cur_filled - 1]] >= 
						cur_alpha[cur_order[k]] - zR * cur_beta[cur_order[k]] - FractionalSolverEpsilon)
						break;
					if (cur_alpha[cur_order[cur_filled - 1]] - cur_range[cur_filled - 1] * cur_beta[cur_order[cur_filled - 1]] - FractionalSolverEpsilon <= 
						cur_alpha[cur_order[k]] - cur_range[cur_filled - 1] * cur_beta[cur_order[k]])
					{
						cur_filled --;
						continue;
					}
					cur_range[cur_filled] = (cur_alpha[cur_order[k]] - cur_alpha[cur_order[cur_filled - 1]]) / (cur_beta[cur_order[k]] - cur_beta[cur_order[cur_filled - 1]]);
					assert(zL <= cur_range[cur_filled] && cur_range[cur_filled] <= zR);
					cur_order[cur_filled ++] = cur_order[k];
					break;
				}
				if (cur_filled == 0)
				{
					cur_range[0] = zL;
					cur_order[cur_filled ++] = cur_order[k];
				}
			}
			for (int k = 1; k < cur_filled; k ++)
				break_point[nbreak_point ++] = cur_range[k];
		}

		break_point[nbreak_point ++] = zL;
		break_point[nbreak_point ++] = zR;
		assert(nbreak_point < n_max * (C + 1) + 2);
		std::sort(break_point, break_point + nbreak_point);
		nbreak_point = (int) (::std::unique(break_point, break_point + nbreak_point) - break_point);

		printf("%d: %d, gap = %e\n", i, nbreak_point, zR - zL);
		if (nbreak_point > 2)
		{
			int index_l = 0, index_r = nbreak_point - 2, index;
			while (index_l <= index_r)
			{
				int index_mid = (index_l + index_r) / 2;
				double v_mid = DPOracle(m, n, C, A, B, a0, b0, break_point[index_mid], dp, NULL);
				if (v_mid >= 0)
				{
					index = index_mid;
					index_l = index_mid + 1;
				}
				else
					index_r = index_mid - 1;
			}
			zL = break_point[index];
			zR = break_point[index + 1];
		}

		double zmid = (zL + zR) * 0.5;
		for (int j = 0; j <= C; j ++)
		{
			int index_max = -1;
			double v_max;
			for (int d = 0; d <= j && d <= n[i]; d ++)
			{
				double v_cur = (last_a[j - d] + A[i][d]) - (last_b[j - d] + B[i][d]) * zmid;
				if (index_max == -1 || v_max < v_cur)
				{
					index_max = d;
					v_max = v_cur; 
				}
			}
			assert(index_max != -1);
			cur_a[j] = last_a[j - index_max] + A[i][index_max];
			cur_b[j] = last_b[j - index_max] + B[i][index_max];
		}
		::std::swap(cur_a, last_a);
		::std::swap(cur_b, last_b);
	}

	z = last_a[C] / last_b[C];
	
	double tmp;
	if (q)
	{
		tmp = DPOracle(m, n, C, A, B, a0, b0, z, dp, prev);
		int p = C;
		for (int i = m; i >= 1; i --)
		{
			q[i - 1] = prev[i][p];
			p -= prev[i][p];
		}
	}
	else
		tmp = DPOracle(m, n, C, A, B, a0, b0, z, dp, NULL);

	if (tmp < -FractionalSolverOptimalityEpsilon || tmp > FractionalSolverOptimalityEpsilon)
		fprintf(stderr, "[FractionalSolver] WARNING: The result may be inaccurate due to precision issues.\n");

	delete last_a;
	delete last_b;
	delete cur_a;
	delete cur_b;
	for (int i = 0; i <= m; i ++)
	{
		delete dp[i];
		delete prev[i];
	}
	delete dp;
	delete prev;
	delete break_point;
	delete cur_alpha;
	delete cur_beta;
	delete cur_order;
	delete cur_range;
}

#ifndef _USE_CPLEX
void MaximizeFraction_CPLEX(int m, int* n, int C, double** A, double** B, double a0, double b0, double& z, int* q)
{
	z = -1;
}
#else
void MaximizeFraction_CPLEX(int m, int* n, int C, double** A, double** B, double a0, double b0, double& z, int* q)
{
	int totElement = 0;
	int n_max = 0;
	for (int i = 0; i < m; i ++)
	{
		totElement += n[i];
		n_max = std::max(n_max, n[i]);
	}
	C = std::max(C, 0);
	C = std::max(C, totElement);

	IloEnv env;
	try
	{
		IloModel model(env);
		IloNumVar _z(env, -IloInfinity, IloInfinity);
		IloNumVarArray theta(env, (m + 1) * (C + 1), -IloInfinity, IloInfinity);
		
		model.add(IloMinimize(env, _z));
		for (int j = 0; j <= C; j ++)
			model.add(theta[j] == a0 - b0 * _z);
		for (int i = 1; i <= m; i ++)
			for (int j = 0; j <= C; j ++)
				for (int d = 0; d <= j && d <= n[i]; d ++)
					model.add(theta[i * (C + 1) + j] >= theta[(i - 1) * (C + 1) + (j - d)] + A[i - 1][d] - B[i - 1][d] * _z);
		model.add(theta[m * (C + 1) + C] <= 0);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		if (cplex.solve() == IloFalse)
			printf("CPLEX Error: No Solution?\n");
		else
			z = cplex.getValue(_z);
		env.end();
	}
	catch (IloException& ex)
	{
		::std::cerr << "Error: " << ex << ::std::endl;
		env.end();
	}
	catch (...)
	{
		::std::cerr << "Error" << ::std::endl;
		env.end();
	}
}
#endif
