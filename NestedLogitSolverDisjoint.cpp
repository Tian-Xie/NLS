/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitSolver.h"

#include <algorithm>
#include <utility>
#include <math.h>
#include <time.h>

const double AssortmentSolverEpsilon = 1e-8;
const double AssortmentSolverInfinity = 1e100;

struct TDelta
{
	double u, Intercept, Slope;
	TDelta() {}
	TDelta(double _u, double _Intercept, double _Slope) : u(_u), Intercept(_Intercept), Slope(_Slope) {}
	bool operator < (const TDelta& b)
	{
		return u < b.u;
	}
};

double NestedLogitSolver::SolveDisjoint(vector <int> C, bool renew_piece)
{
	int m = Model.nNest;
	for (int i = 0; i < m; i ++)
	{
		C[i] = std::max(C[i], 0);
		C[i] = std::min(C[i], Model.Nest[i].nProduct);
	}
	if (renew_piece)
	{
		for (int i = 0; i < m; i ++)
			Model.Nest[i].GeneratePiecewiseLinear(C[i]);
	}

	// Merge piecewise-linear functions
	int nDelta = 0;
	for (int i = 0; i < m; i ++)
		nDelta += Model.Nest[i].nPiece[C[i]] - 1;

	TDelta* Delta = new TDelta[nDelta];
	double Cur_Intercept = 0, Cur_Slope = Model.V0;
	nDelta = 0;
	for (int i = 0; i < m; i ++)
	{
		TCandidate* lines = Model.Nest[i].Candidate[C[i]];
		int nlines = Model.Nest[i].nCandidate[C[i]];
		TPiece* Piece = Model.Nest[i].Piece[C[i]];
		int nPiece = Model.Nest[i].nPiece[C[i]];

		// Make deltas
		Cur_Intercept += lines[Piece[0].first].Intercept;
		Cur_Slope += lines[Piece[0].first].Slope;
		for (int j = 1; j < nPiece; j ++)
		{
			Delta[nDelta ++] = TDelta(
				Piece[j].second, 
				lines[Piece[j].first].Intercept - lines[Piece[j - 1].first].Intercept, 
				lines[Piece[j].first].Slope - lines[Piece[j - 1].first].Slope
			);
		}
	}

	// Sort and scan
	std::sort(Delta, Delta + nDelta);
	double Objective = -AssortmentSolverInfinity;
	if (fabs(Cur_Slope) > AssortmentSolverEpsilon) Objective = std::max(Objective, Cur_Intercept / Cur_Slope);
	for (int d = 0; d < nDelta; d ++)
	{
		Cur_Intercept += Delta[d].Intercept;
		Cur_Slope += Delta[d].Slope;
		if (fabs(Cur_Slope) > AssortmentSolverEpsilon) Objective = std::max(Objective, Cur_Intercept / Cur_Slope);
	}

	delete Delta;
	OptimalZ = Objective;
	OptimalCon = C;
	return OptimalZ;
}
