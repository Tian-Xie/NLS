/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "NestedLogitModel.h"
#include <random>
#include <time.h>

bool TNest::InsertProduct(TProduct p)
{
	if (p.R < 0 || p.V < 0) return false;
	nProduct ++;
	Product.push_back(p);
	return true;
}

void TNest::Clear()
{
	nProduct = 0;
	Product.clear();
	if (Candidate)
	{
		for (int nnc = 0; nnc < nnCandidate; nnc ++)
			delete Candidate[nnc];
		delete Candidate;
		delete nCandidate;
		nnCandidate = 0;
		nCandidate = NULL;
		Candidate = NULL;
	}
	if (Piece)
	{
		for (int nnp = 0; nnp < nnPiece; nnp ++)
			if (nPiece[nnp]) delete Piece[nnp];
		delete Piece;
		delete nPiece;
		nnPiece = 0;
		nPiece = NULL;
		Piece = NULL;
	}
}

TNest::~TNest()
{
	Clear();
}

NestedLogitModel::NestedLogitModel()
{
	Reset();
	SetV0(1.0);
}

NestedLogitModel::NestedLogitModel(double _V0)
{
	Reset();
	SetV0(_V0);
}

NestedLogitModel::NestedLogitModel(const NestedLogitModel& MOD)
{
	nNest = MOD.nNest;
	Nest = MOD.Nest;
	V0 = MOD.V0;
}

void NestedLogitModel::SetV0(double _V0)
{
	V0 = _V0;
}

void NestedLogitModel::Reset()
{
	for (int t = 0; t < Nest.size(); t ++)
		Nest[t].Clear();
	nNest = 0;
	Nest.clear();
}

int NestedLogitModel::InsertNest(double _Gamma)
{
	if (_Gamma <= 0 || _Gamma > 1) return -1;
	Nest.push_back(TNest(_Gamma));
	return nNest ++;
}

int NestedLogitModel::InsertProduct(int NestID, double v, double r)
{
	if (NestID >= nNest) return -1;
	if (v < 0) return -1;
	if (r < 0) return -1;
	Nest[NestID].InsertProduct(TProduct(v, r));
	return Nest[NestID].nProduct - 1;
}

void NestedLogitModel::RandomCase(int nNest, int nProduct, time_t seed)
{
	Reset();
	SetV0(nNest * 10.0);
	std::uniform_real_distribution <double> RandGamma(0.25, 0.75);
	std::uniform_real_distribution <double> RandC(0.0, 1.0);
	std::uniform_real_distribution <double> RandXY(0.75, 1.25);
	
	std::mt19937 Gen(seed);
	for (int i = 0; i < nNest; i ++)
	{
		InsertNest(RandGamma(Gen));
		for (int j = 0; j < nProduct; j ++)
		{
			double C = RandC(Gen);
			double r = 10.0 * C * C * RandXY(Gen);
			double v = 10.0 * (1.0 - C) * RandXY(Gen);
			InsertProduct(i, v, r);
		}
	}
}

NestedLogitModel::~NestedLogitModel()
{
	Reset();
}

