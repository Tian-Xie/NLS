/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include "FenwickTree.h"

FenwickTree::FenwickTree(int _n)
{
	n = _n;
	S = new double[n + 1];
	for (int i = 0; i <= n; i ++)
		S[i] = 0;
}

void FenwickTree::MakeChange(int index, double Delta)
{
	for (index ++; index <= n; index += (index & (-index)))
		S[index] += Delta;
}

double FenwickTree::GetPrefixSum(int index)
{
	double Ret = 0;
	for (index ++; index > 0; index -= (index & (-index)))
		Ret += S[index];
	return Ret;
}

FenwickTree::~FenwickTree()
{
	delete S;
}
