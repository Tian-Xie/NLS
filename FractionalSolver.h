/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

void MaximizeFraction(int m, int* n, int C, double** A, double** B, double a0, double b0, double& z, int* q);
void MaximizeFraction_CPLEX(int m, int* n, int C, double** A, double** B, double a0, double b0, double& z, int* q);
