#include "NestedLogitSolver.h"
#include <vector>

using std::vector;

void TestDisjoint(FILE* out, int m, int n, vector <double> _ratio, int useCPLEX);
void TestJoint(FILE* out, int m, int n, vector <double> _ratio, int useCPLEX);
