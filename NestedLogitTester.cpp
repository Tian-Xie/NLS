#include "NestedLogitTester.h"
#include <time.h>
#include <math.h>
#include <vector>

using std::vector;

void TestDisjoint(FILE* out, int m, int n, vector <double> _ratio, int useCPLEX)
{
	NestedLogitModel Model;
	int seed = time(0);
	Model.RandomCase(m, n, seed);

	clock_t startTime, endTime;
	double sec, ans;

	if (out) fprintf(out, "------ Disjoint: m = %d, n = %d, seed = %d, ------\n", m, n, seed);
	fprintf(stderr, "------ Disjoint: m = %d, n = %d, seed = %d ------\n", m, n, seed);

	startTime = clock();
	Model.GenerateAllCandidateSet();
	endTime = clock();
	sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;

	if (out) fprintf(out, "GenerateAllCandidateSet: %.4lf\n", sec);
	fprintf(stderr, "GenerateAllCandidateSet: %.4lf\n", sec);
	
	NestedLogitSolver Solver(Model);

	for (auto ratio: _ratio)
	{
		if (out) fprintf(out, "------ ratio: %.2lf ------\n", ratio);
		fprintf(stderr, "------ ratio: %.2lf ------\n", ratio);

		vector <int> C(m, (int) ceil(n * ratio));
		startTime = clock();
		ans = Solver.SolveDisjoint(C);
		endTime = clock();
		sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;
		if (out) fprintf(out, "%.10lf \t SolveDisjoint: %.4lf\n", ans, sec);
		fprintf(stderr, "%.10lf \t SolveDisjoint: %.4lf\n", ans, sec);

		if (useCPLEX)
		{
			startTime = clock();
			ans = Solver.SolveDisjointCPLEX(C);
			endTime = clock();
			sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;
			if (out) fprintf(out, "%.10lf \t SolveDisjointCPLEX: %.4lf\n", ans, sec);
			fprintf(stderr, "%.10lf \t SolveDisjointCPLEX: %.4lf\n", ans, sec);
		}
	}

	if (out) fprintf(out, "------ Finish ------\n");
	fprintf(stderr, "------ Finish ------\n");
}

void TestJoint(FILE* out, int m, int n, vector <double> _ratio, int useCPLEX)
{
	NestedLogitModel Model;
	int seed = time(0);
	Model.RandomCase(m, n, seed);

	clock_t startTime, endTime;
	double sec, ans;

	if (out) fprintf(out, "------ Joint: m = %d, n = %d, seed = %d ------\n", m, n, seed);
	fprintf(stderr, "------ Joint: m = %d, n = %d, seed = %d ------\n", m, n, seed);

	startTime = clock();
	Model.GenerateAllCandidateSet();
	endTime = clock();
	sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;

	if (out) fprintf(out, "GenerateAllCandidateSet: %.4lf\n", sec);
	fprintf(stderr, "GenerateAllCandidateSet: %.4lf\n", sec);
	
	NestedLogitSolver Solver(Model);

	for (auto ratio: _ratio)
	{
		if (out) fprintf(out, "------ ratio: %.2lf ------\n", ratio);
		fprintf(stderr, "------ ratio: %.2lf ------\n", ratio);
	
		int C = (int) ceil(m * n * ratio);
	
		startTime = clock();
		ans = Solver.SolveJoint(C);
		endTime = clock();
		sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;
		if (out) fprintf(out, "%.10lf \t SolveJoint: %.4lf\n", ans, sec);
		fprintf(stderr, "%.10lf \t SolveJoint: %.4lf\n", ans, sec);

		if (useCPLEX)
		{
			startTime = clock();
			ans = Solver.SolveJointCPLEX(C);
			endTime = clock();
			sec = (double) (endTime - startTime) / CLOCKS_PER_SEC;
			if (out) fprintf(out, "%.10lf \t SolveJointCPLEX: %.4lf\n", ans, sec);
			fprintf(stderr, "%.10lf \t SolveJointCPLEX: %.4lf\n", ans, sec);
		}
	}

	if (out) fprintf(out, "------ Finish ------\n");
	fprintf(stderr, "------ Finish ------\n");
}

