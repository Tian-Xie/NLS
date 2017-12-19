/*
	NLS: NestedLogitSolver
	Constrained Assortment Optimization Solver for Nested Logit Model
	Author: Tian Xie (SHUFE)
	Date: Oct 19, 2017
*/

#include <vector>
#include <utility>
#include <stdio.h>
#include <time.h>
#include "NestedLogitModel.h"
#include "NestedLogitSolver.h"
#include "NestedLogitTester.h"

using namespace std;

int main()
{
	FILE* out = fopen("output.txt", "w");

	try
	{
		vector <double> ratio;
		ratio.push_back(0.1);
		ratio.push_back(0.15);
		ratio.push_back(0.2);

		int nRM = 1;
		int RM[] = {1000};
		int nRN = 1;
		int RN[] = {100};

		vector < pair <int, int> > Case;
		for (int m = 0; m < nRM; m ++)
			for (int n = 0; n < nRN; n ++)
				Case.push_back(make_pair(RM[m], RN[n]));
		
		for (auto c: Case)
		{
			TestJoint(out, c.first, c.second, ratio, false);
			fflush(out);
		}
	}
	catch (const bad_alloc& e)
	{
		fprintf(stderr, "bad_alloc!\n");
	}
	
	
	auto ans = Solver.GetOptimalAssortment();
	for (int i = 0; i < ans.size(); i ++)
	{
		printf("Nest #%d: [%d chosen] [ ", i, (int) ans[i].size());
		for (auto p: ans[i]) printf("%d ", p);
		printf("]\n");
	}
	printf("\n");
	

	return 0;
}
