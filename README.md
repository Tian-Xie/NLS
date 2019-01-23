# NLS: NestedLogitSolver

Constrained Assortment Optimization Solver for Nested Logit Model


### Author

Tian Xie (Research Center for Management Science and Information Analytics, Shanghai University of Finance and Economics)

### What is this?

This code is an implementation of a fast algorithm which solves constrained assortment optimization problem (Gallego and Topaloglu, 2014; Feldman and Topaloglu, 2015) without linear programming. 

If you use this code, please cite: 

Tian Xie, and Dongdong Ge. 2018. A Tractable Discrete Fractional Programming: Application to Constrained Assortment Optimization. Journal of Combinatorial Optimization 36(2): 400-415.

### How to use it?

Step 1. The user should manually create an instance of `NestedLogitModel` and set related model information. The function `NestedLogitModel::RandomCase` in `NestedLogitModel.cpp` gives a detailed information about how to initialize, create nests and products. 

Step 2. Call method `GenerateAllCandidateSet` (or equivalent operators) for the model instance to generate candidate sets. The algorithm is stated in Appendix A of our paper. A code file named `NestedLogitTester.cpp` demonstrates how to do this.

Step 3. Create an instance of `NestedLogitSolver` and bind your `NestedLogitModel` instance to it. Then use `NestedLogitSolver::SolveDisjoint` or `NestedLogitSolver::SolveJoint` to obtain optimal expected revenue. `NestedLogitTester.cpp` is a reference. If you want an optimal assortment, you can see how to use `NestedLogitSolver::GetOptimalAssortment` in `Demo.cpp`. 

Note that an independent FractionalSolver can be used to solve the fractional programming problem seperately. 

### External Dependencies

The problems are also applicable for directly solving linear programming, so we also prepare a CPLEX interface in `NestedLogitCPLEXInterface.cpp` and you can see how to call it in `NestedLogitTester.cpp`. If you want to have a comparison, please link your CPLEX solver yourself and add 
`#define _USE_CPLEX`
to this code before your compilation to activate the functions in `NestedLogitCPLEXInterface.cpp`. 

### References

[1] Gallego G., and Topaloglu H. 2014. Constrained assortment optimization for the nested logit model. Management Science 60(10): 2583-2601. 

[2] Feldman J. B., and Topaloglu H. 2015. Capacity Constraints Across Nests in Assortment Optimization Under the Nested Logit Model. Operations Research 63(4): 812-822.
