/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#ifndef _FENWICKTREE_H
#define _FENWICKTREE_H

struct FenwickTree
{
	int n;
	double* S;

	FenwickTree(int _n);
	void MakeChange(int index, double Delta);
	double GetPrefixSum(int index);
	~FenwickTree();
};

#endif
