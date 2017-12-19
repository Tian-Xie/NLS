/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#ifndef _NESTEDLOGITSOLVER_H
#define _NESTEDLOGITSOLVER_H

#include "NestedLogitModel.h"

struct NestedLogitSolver
{
	NestedLogitModel& Model;

	vector <int> OptimalCon;
	double OptimalZ;

	NestedLogitSolver(NestedLogitModel& _Model) : Model(_Model)
	{
		OptimalZ = -1;
		OptimalCon.clear();
	};

	// Private
	double JointDPOracle(double z, int C, double* working_tmp);

	double SolveDisjoint(vector <int> C, bool renew_piece = true);
	double SolveJoint(int C, bool renew_piece = true);
	vector < vector <int> > GetOptimalAssortment();

	// For double-check only
	double SolveDisjointCPLEX(vector <int> C);
	double SolveJointCPLEX(int C);
};

#endif
