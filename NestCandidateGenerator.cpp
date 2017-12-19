/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitModel.h"

#include <utility>
#include <algorithm>
#include <assert.h>
#include "FenwickTree.h"

const double CandidateGeneratorEpsilon = 1e-8;

struct TCrossPoint
{
	int j1, j2;
	TCrossPoint() {}
	TCrossPoint(int _j1, int _j2) : j1(_j1), j2(_j2) {}
};

bool TCandidate::operator < (const TCandidate& b)
{
	if (fabs(Slope - b.Slope) > CandidateGeneratorEpsilon) return Slope > b.Slope;
	return Intercept > b.Intercept;
}

TCandidate MakeCandidate(double V, double RV, double gamma)
{
	if (V < CandidateGeneratorEpsilon)
		return TCandidate(V, RV, 0, 0);
	else
		return TCandidate(
			V, 
			RV, 
			pow(V, gamma - 1.0) * RV, 
			pow(V, gamma)
		);
}

void TNest::GenerateCandidateSet()
{
	int n = nProduct;

	std::pair <double, double>* lines = new std::pair <double, double> [n + 1];
	lines[0] = std::make_pair(0, 0); // auxiliary line 0
	for (int j = 0; j < n; j ++)
		lines[j + 1] = std::make_pair(Product[j].V, Product[j].V * Product[j].R);
	std::sort(lines, lines + (n + 1));

	// Generate all crosspoints
	TCrossPoint* Cross = new TCrossPoint[(n + 1) * n / 2];
	int nCross = 0;
	for (int j1 = 0; j1 <= n; j1 ++)
		for (int j2 = j1 + 1; j2 <= n; j2 ++)
		{
			// Ignore parallel lines
			if (fabs(lines[j1].first - lines[j2].first) < CandidateGeneratorEpsilon) continue;
			//double u = (lines[j1].second - lines[j2].second) / (lines[j1].first - lines[j2].first);
			Cross[nCross ++] = TCrossPoint(j1, j2);
		}
	std::sort(Cross, Cross + nCross, [&](const TCrossPoint& a, const TCrossPoint& b) -> bool
	{
		double v = (lines[a.j2].second - lines[a.j1].second) * (lines[b.j2].first - lines[b.j1].first) - 
			(lines[a.j2].first - lines[a.j1].first) * (lines[b.j2].second - lines[b.j1].second);
		if (fabs(v) > CandidateGeneratorEpsilon) return v < 0;
		if (a.j1 != b.j1) return a.j1 > b.j1;
		return a.j2 > b.j2;
	});

	// Maintain the prefix sums by Fenwick Trees
	FenwickTree TreeV(n + 1), TreeRV(n + 1);
	for (int j = 0; j <= n; j ++)
	{
		TreeV.MakeChange(j, lines[n - j].first);
		TreeRV.MakeChange(j, lines[n - j].second);
	}

	// Initial candidate sets
	if (Candidate)
	{
		for (int nnc = 0; nnc < nnCandidate; nnc ++)
			delete Candidate[nnc];
		delete Candidate;
		delete nCandidate;
		Candidate = NULL;
	}
	nnCandidate = n + 1;
	nCandidate = new int[nnCandidate];
	Candidate = new TCandidate*[nnCandidate];

	int* index = new int[n + 1];
	// Calculate the size of candidate set
	for (int j = 0; j <= n; j ++)
	{
		index[j] = n - j;
		nCandidate[j] = 1;
	}
	for (int cp_index = 0; cp_index < nCross; cp_index ++)
	{
		TCrossPoint& cp = Cross[cp_index];
		if (cp.j1 == 0)
			for (int c = index[0]; c <= n; c ++) nCandidate[c] ++;
		else if (index[cp.j1] < index[0])
			nCandidate[index[cp.j1]] ++;
		std::swap(index[cp.j1], index[cp.j2]);
	}
	for (int j = 0; j <= n; j ++)
		Candidate[j] = new TCandidate[nCandidate[j]];

	// Reset
	for (int j = 0; j <= n; j ++)
	{
		index[j] = n - j;
		nCandidate[j] = 1;
	}
	for (int j = 0; j <= n; j ++)
	{
		Candidate[j][0] = MakeCandidate(
			TreeV.GetPrefixSum(std::min(j - 1, index[0])), 
			TreeRV.GetPrefixSum(std::min(j - 1, index[0])), 
			Gamma
		);
	}

	// Iterate over all crosspoints
	for (int cp_index = 0; cp_index < nCross; cp_index ++)
	{
		TCrossPoint& cp = Cross[cp_index];
		// j1 > j2 <=> j1 < j2
		assert(index[cp.j2] + 1 == index[cp.j1]);
		
		TreeV.MakeChange(index[cp.j1], -lines[cp.j1].first);
		TreeV.MakeChange(index[cp.j1], lines[cp.j2].first);
		TreeV.MakeChange(index[cp.j2], -lines[cp.j2].first);
		TreeV.MakeChange(index[cp.j2], lines[cp.j1].first);
		TreeRV.MakeChange(index[cp.j1], -lines[cp.j1].second);
		TreeRV.MakeChange(index[cp.j1], lines[cp.j2].second);
		TreeRV.MakeChange(index[cp.j2], -lines[cp.j2].second);
		TreeRV.MakeChange(index[cp.j2], lines[cp.j1].second);

		if (cp.j1 == 0)
		{
			double tmpV = TreeV.GetPrefixSum(index[0] - 1);
			double tmpRV = TreeRV.GetPrefixSum(index[0] - 1);
			for (int c = index[0]; c <= n; c ++)
				Candidate[c][nCandidate[c] ++] = MakeCandidate(tmpV, tmpRV, Gamma);
		}
		else if (index[cp.j1] < index[0])
		{
			Candidate[index[cp.j1]][nCandidate[index[cp.j1]] ++] = MakeCandidate(
				TreeV.GetPrefixSum(index[cp.j1] - 1), 
				TreeRV.GetPrefixSum(index[cp.j1] - 1), 
				Gamma
			);
		}
		std::swap(index[cp.j1], index[cp.j2]);
	}

	for (int j = 0; j <= n; j ++)
		std::sort(Candidate[j], Candidate[j] + nCandidate[j]);

	delete index;
	delete lines;
	delete Cross;
}

void NestedLogitModel::GenerateAllCandidateSet()
{
	for (int i = 0; i < nNest; i ++)
		Nest[i].GenerateCandidateSet();
}
