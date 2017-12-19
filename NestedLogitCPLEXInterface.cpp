/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitSolver.h"

#ifdef _USE_CPLEX

#include <iostream>
#include <math.h>
#include <ilcplex/ilocplex.h>

double NestedLogitSolver::SolveDisjointCPLEX(vector <int> C)
{
	int m = Model.nNest;
	for (int i = 0; i < m; i ++)
	{
		C[i] = std::max(C[i], 0);
		C[i] = std::min(C[i], Model.Nest[i].nProduct);
	}

	IloEnv env;
	try
	{
		IloModel model(env);
		IloNumVarArray Vy(env, Model.nNest, -IloInfinity, IloInfinity);
		IloNumVar Vz(env, -IloInfinity, IloInfinity);
		
		model.add(IloMinimize(env, Vz));
		model.add(Model.V0 * Vz >= IloSum(Vy));
		for (int i = 0; i < m; i ++)
		{
			int nc = Model.Nest[i].nCandidate[C[i]];
			TCandidate* pc = Model.Nest[i].Candidate[C[i]];
			for (int j = 0; j < nc; j ++)
			{
				TCandidate& can = pc[j];
				double canR = (can.V == 0) ? 0 : (can.RV / can.V);
				model.add(Vy[i] >= pow(can.V, Model.Nest[i].Gamma) * (canR - Vz));
			}
		}
		
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::TiLim, 7200.0);
		cplex.solve();
		OptimalZ = cplex.getObjValue();
		env.end();
	}
	catch (IloException& ex)
	{
		std::cerr << "Error: " << ex << std::endl;
		OptimalZ = -1.0;
		env.end();
	}
	catch (...)
	{
		std::cerr << "Unknown Error" << std::endl;
		OptimalZ = -1.0;
		env.end();
	}
	return OptimalZ;
}

#define IND(i, q) ((i + 1) * (C + 1) + (q))

double NestedLogitSolver::SolveJointCPLEX(int C)
{
	int m = Model.nNest;
	IloEnv env;
	try
	{
		IloModel model(env);
		IloNumVarArray Vtheta(env, (m + 1) * (C + 1), -IloInfinity, IloInfinity);
		IloNumVar Vz(env, -IloInfinity, IloInfinity);
		
		model.add(IloMinimize(env, Vz));
		model.add(Vtheta[IND(m - 1, C)] <= 0);
		for (int i = 0; i < m; i ++)
		{
			int n = Model.Nest[i].nProduct;
			for (int q = 0; q <= C; q ++)
			{
				for (int d = 0; d <= n && d <= q; d ++)
				{
					int nc = Model.Nest[i].nCandidate[d];
					TCandidate* pc = Model.Nest[i].Candidate[d];
					for (int j = 0; j < nc; j ++)
					{
						TCandidate& can = pc[j];
						model.add(Vtheta[IND(i, q)] >= Vtheta[IND(i - 1, q - d)] + can.Intercept - can.Slope * Vz);
					}
				}
			}
		}
		for (int q = 0; q <= C; q ++)
			model.add(Vtheta[IND(-1, q)] == -Model.V0 * Vz);
		
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::TiLim, 7200.0);
		cplex.solve();
		OptimalZ = cplex.getObjValue();
		env.end();
	}
	catch (IloException& ex)
	{
		std::cerr << "Error: " << ex << std::endl;
		OptimalZ = -1.0;
		env.end();
	}
	catch (...)
	{
		std::cerr << "Unknown Error" << std::endl;
		OptimalZ = -1.0;
		env.end();
	}
	return OptimalZ;
}

#else

double NestedLogitSolver::SolveDisjointCPLEX(vector <int> C)
{
	return -1;
}

double NestedLogitSolver::SolveJointCPLEX(int C)
{
	return -1;
}

#endif
