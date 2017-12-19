/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitSolver.h"

#include <algorithm>
#include <assert.h>
#include "FractionalSolver.h"

const double AssortmentSolverEpsilon = 1e-8;
const double AssortmentSolverInfinity = 1e100;

double NestedLogitSolver::JointDPOracle(double z, int C, double* working_tmp)
{
	double* Prev = working_tmp;
	double* Cur = working_tmp + (C + 1);
	double* Val = working_tmp + (2 * (C + 1));
	
	int m = Model.nNest;
	for (int c = 0; c <= C; c ++)
		Cur[c] = -Model.V0 * z;

	for (int i = 0; i < m; i ++)
	{
		std::swap(Cur, Prev);
		TNest& curnest = Model.Nest[i];
		int nprod = curnest.nProduct;
		for (int j = 0; j <= nprod; j ++)
		{
			double tmp = -AssortmentSolverInfinity;
			int n_can_arr = curnest.nCandidate[j];
			TCandidate* can_arr = curnest.Candidate[j];
			for (int k = 0; k < n_can_arr; k ++)
				if (can_arr[k].Intercept - can_arr[k].Slope * z > tmp)
					tmp = can_arr[k].Intercept - can_arr[k].Slope * z;
			Val[j] = tmp;
		}
		for (int c = 0; c <= C; c ++)
		{
			double tmp = -AssortmentSolverInfinity;
			for (int j = 0; j <= nprod && j <= c; j ++)
				if (Prev[c - j] + Val[j] > tmp)
					tmp = Prev[c - j] + Val[j];
			Cur[c] = tmp;
		}
	}
	return Cur[C];
}

double NestedLogitSolver::SolveJoint(int C, bool renew_piece)
{
	int m = Model.nNest;
	int totProduct = 0;
	int n_max = 0;
	int* n = new int[m];
	for (int i = 0; i < m; i ++)
	{
		n[i] = Model.Nest[i].nProduct;
		totProduct += n[i];
		n_max = std::max(n_max, n[i]);
	}
	C = std::max(C, 0);
	C = std::min(C, totProduct);
	for (int i = 0; i < m; i ++)
		n[i] = std::min(n[i], C);

	if (renew_piece)
	{
		for (int i = 0; i < m; i ++)
			for (int j = 0; j <= n[i]; j ++)
				Model.Nest[i].GeneratePiecewiseLinear(j);
	}

	// Collect all breakpoints
	int n_break_points = 0;
	for (int i = 0; i < m; i ++)
		for (int j = 0; j <= n[i]; j ++)
			n_break_points += Model.Nest[i].nPiece[j];
	double* break_points = new double[n_break_points + 1];
	n_break_points = 0;
	for (int i = 0; i < m; i ++)
	{
		TNest& curnest = Model.Nest[i];
		for (int j = 0; j <= n[i]; j ++)
			for (int k = 0; k < curnest.nPiece[j]; k ++)
				break_points[n_break_points ++] = curnest.Piece[j][k].second;
	}
	break_points[n_break_points ++] = AssortmentSolverInfinity;
	std::sort(break_points, break_points + n_break_points);

	// Binary search
	int L = 0, R = n_break_points - 2, Ptr = -1;
	// Temporary
	double* working_tmp = new double[2 * (C + 1) + n_max + 1];

	while (L <= R)
	{
		int Mid = (L + R) / 2;
		double obj = JointDPOracle(break_points[Mid], C, working_tmp);
		if (obj >= 0)
		{
			Ptr = Mid;
			L = Mid + 1;
		}
		else
			R = Mid - 1;
	}
	assert(Ptr != -1);
	double midz = (break_points[Ptr] + break_points[Ptr + 1]) * 0.5;

	// Construct an instance for FractionalSolver
	double** FractionalSolver_A = new double*[m];
	double** FractionalSolver_B = new double*[m];
	for (int i = 0; i < m; i ++)
	{
		FractionalSolver_A[i] = new double[n[i] + 1];
		FractionalSolver_B[i] = new double[n[i] + 1];
		TNest& curnest = Model.Nest[i];
		for (int j = 0; j <= n[i]; j ++)
		{
			TCandidate* can = curnest.Candidate[j];
			double curbest = -AssortmentSolverInfinity;
			int curbest_ptr = -1;
			for (int k = 0; k < curnest.nCandidate[j]; k ++)
				if (can[k].Intercept - can[k].Slope * midz > curbest)
				{
					curbest = can[k].Intercept - can[k].Slope * midz;
					curbest_ptr = k;
				}
			FractionalSolver_A[i][j] = can[curbest_ptr].Intercept;
			FractionalSolver_B[i][j] = can[curbest_ptr].Slope;
		}
	}
	// Solve by FractionalSolver
	double FractionalSolver_z = -1;
	int* FractionalSolver_q = new int[m];
	MaximizeFraction(m, n, C, FractionalSolver_A, FractionalSolver_B, 0, Model.V0, FractionalSolver_z, FractionalSolver_q);
	// Refine solution by disjoint solver
	OptimalCon.clear();
	OptimalCon.resize(m);
	for (int i = 0; i < m; i ++)
		OptimalCon[i] = FractionalSolver_q[i];
	OptimalZ = -1;
	SolveDisjoint(OptimalCon, false);
	// Precision Error Detection
	assert(fabs(FractionalSolver_z - OptimalZ) < AssortmentSolverEpsilon);

	delete working_tmp;
	delete break_points;
	for (int i = 0; i < m; i ++)
	{
		delete FractionalSolver_A[i];
		delete FractionalSolver_B[i];
	}
	delete FractionalSolver_A;
	delete FractionalSolver_B;
	delete FractionalSolver_q;
	delete n;
	return OptimalZ;
}
