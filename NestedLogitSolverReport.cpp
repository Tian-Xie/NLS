/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitSolver.h"

#include <algorithm>
#include <assert.h>

const double ReporterEpsilon = 1e-8;
const double ReporterInfinity = 1e100;

vector < vector <int> > NestedLogitSolver::GetOptimalAssortment()
{
	vector < vector <int> > Res;
	double z = OptimalZ;
	int m = Model.nNest;
	for (int i = 0; i < m; i ++)
	{
		TNest& curnest = Model.Nest[i];
		int ci = OptimalCon[i];
		int ncan = curnest.nCandidate[ci];
		TCandidate* can = curnest.Candidate[ci];

		double best = -ReporterInfinity;
		int best_index = -1;
		for (int j = 0; j < ncan; j ++)
			if (can[j].Intercept - can[j].Slope * z > best)
			{
				best = can[j].Intercept - can[j].Slope * z;
				best_index = j;
			}
		
		double Rstar = 0;
		if (can[best_index].V > ReporterEpsilon) Rstar = can[best_index].RV / can[best_index].V;
		double u = std::max(z, curnest.Gamma * z + (1.0 - curnest.Gamma) * Rstar);
		
		vector <int> Res_Curnest(curnest.nProduct);
		for (int j = 0; j < curnest.nProduct; j ++)
			Res_Curnest[j] = j;
		vector <TProduct>& P = curnest.Product;
		std::sort(Res_Curnest.begin(), Res_Curnest.end(), [&](const int& a, const int& b) -> bool
		{
			return P[a].V * (P[a].R - u) > P[b].V * (P[b].R - u);
		});
		while (ci > 0 && P[Res_Curnest[ci - 1]].V * (P[Res_Curnest[ci - 1]].R - u) <= ReporterEpsilon)
			ci --;
		Res_Curnest.resize(ci);
		sort(Res_Curnest.begin(), Res_Curnest.end());
		Res.push_back(Res_Curnest);

		double tmpRV = 0, tmpV = 0;
		for (int ind: Res_Curnest)
		{
			tmpRV += P[ind].R * P[ind].V;
			tmpV += P[ind].V;
		}
		assert(fabs(tmpRV - can[best_index].RV) < ReporterEpsilon);
		assert(fabs(tmpV - can[best_index].V) < ReporterEpsilon);
	}
	return Res;
}
