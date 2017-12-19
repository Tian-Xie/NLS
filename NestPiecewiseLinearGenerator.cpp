/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitModel.h"

const double PiecewiseLinearGeneratorEpsilon = 1e-8;
const double PiecewiseLinearGeneratorInfinity = 1e100;

void TNest::GeneratePiecewiseLinear(int C = -1)
{
	int n = nProduct;
	int* new_nPiece = new int[n + 1];
	TPiece** new_Piece = new TPiece*[n + 1];
	for (int nnp = 0; nnp <= n; nnp ++)
	{
		new_nPiece[nnp] = 0;
		new_Piece[nnp] = NULL;
	}
	if (Piece)
	{
		for (int nnp = 0; nnp < nnPiece && nnp <= n; nnp ++)
			if (nPiece[nnp])
			{
				new_nPiece[nnp] = nPiece[nnp];
				new_Piece[nnp] = Piece[nnp];
			}
		for (int nnp = n + 1; nnp < nnPiece; nnp ++)
			if (nPiece[nnp])
				delete Piece[nnp];
		delete Piece;
		delete nPiece;
		Piece = NULL;
	}
	Piece = new_Piece;
	nPiece = new_nPiece;
	nnPiece = n + 1;

	int cL = 0, cR = nProduct;
	if (C != -1)
	{
		C = std::max(C, 0);
		C = std::min(C, nProduct);
		cL = cR = C;
	}

	int nlines_max = 0;
	for (int c = cL; c <= cR; c ++)
		nlines_max = std::max(nlines_max, nCandidate[c]);
	TPiece* tPiece = new TPiece[nlines_max];

	for (int c = cL; c <= cR; c ++)
	{
		TCandidate* lines = Candidate[c];
		int nlines = nCandidate[c];

		int ntPiece = 1;
		tPiece[0] = std::make_pair(0, -PiecewiseLinearGeneratorInfinity);
		for (int j = 1; j < nlines; j ++)
		{
			if (fabs(lines[j].Slope - lines[j - 1].Slope) < PiecewiseLinearGeneratorEpsilon) continue;
			while (ntPiece > 1)
			{
				int jPrime = tPiece[ntPiece - 1].first;
				double x = tPiece[ntPiece - 1].second;
				double val1 = lines[jPrime].Intercept - x * lines[jPrime].Slope;
				double val2 = lines[j].Intercept - x * lines[j].Slope;
				if (val1 <= val2 + PiecewiseLinearGeneratorEpsilon)
					ntPiece --;
				else
					break;
			}
			int jPrime = tPiece[ntPiece - 1].first;
			tPiece[ntPiece ++] = std::make_pair(
				j, 
				(lines[jPrime].Intercept - lines[j].Intercept) / (lines[jPrime].Slope - lines[j].Slope)
			);
		}
		nPiece[c] = ntPiece;
		Piece[c] = new TPiece[ntPiece];
		for (int j = 0; j < ntPiece; j ++)
			Piece[c][j] = tPiece[j];
	}

	delete tPiece;
}

void NestedLogitModel::GenerateAllPiecewiseLinear()
{
	for (int i = 0; i < nNest; i ++)
		Nest[i].GeneratePiecewiseLinear();
}
