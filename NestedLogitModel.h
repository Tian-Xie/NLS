/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#ifndef _NESTEDLOGITMODEL_H
#define _NESTEDLOGITMODEL_H

#include <vector>
#include <utility>
#include <math.h>
#include <time.h>

using std::vector;

struct TProduct
{
	double V, R;
	TProduct(double _V, double _R): V(_V), R(_R) {}
};

struct TCandidate
{
	double V, RV;
	double Intercept, Slope; // Intercept - Slope * u
	TCandidate() {}
	TCandidate(double _V, double _RV, double _Intercept, double _Slope) : V(_V), RV(_RV), Intercept(_Intercept), Slope(_Slope) {}
	bool operator < (const TCandidate& b);
};

typedef std::pair <int, double> TPiece;

struct TNest
{
	int nProduct;
	double Gamma;
	vector <TProduct> Product;
	
	TCandidate** Candidate; // It may cost O(nProduct^2) to storage
	int* nCandidate;
	int nnCandidate;
	TPiece** Piece; // Save the piecewise-linear form
	int* nPiece;
	int nnPiece;

	TNest(double _Gamma): Gamma(_Gamma)
	{
		nProduct = 0;
		Product.clear();
		Candidate = NULL;
		nCandidate = NULL;
		nnCandidate = 0;
		Piece = NULL;
		nPiece = NULL;
		nnPiece = 0;
	}

	bool InsertProduct(TProduct p);
	void Clear();
	void GenerateCandidateSet();
	void GeneratePiecewiseLinear(int C);
	~TNest();
};

struct NestedLogitModel
{
	int nNest; // How many nests?
	vector <TNest> Nest;
	double V0;
	
	NestedLogitModel();
	NestedLogitModel(double _V0);
	NestedLogitModel(const NestedLogitModel& MOD);
	void SetV0(double _V0);
	void Reset();
	int InsertNest(double _Gamma);
	int InsertProduct(int NestID, double v, double r);
	void RandomCase(int nNest, int nProduct, time_t seed);
	void GenerateAllCandidateSet();
	void GenerateAllPiecewiseLinear();
	~NestedLogitModel();
};

#endif
