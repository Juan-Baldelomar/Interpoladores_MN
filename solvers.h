//
// Created by juan on 23/10/20.
//

#ifndef MN_T9_SOLVERS_H
#define MN_T9_SOLVERS_H

#include "math.h"
#include <vector>

using namespace std;
//utilities
void TransposeMatrix(vector<vector<double>>& A, vector<vector<double>>&B);

double RecursiveNestedNewtonForm(vector<double>&a, vector<double>&c, double x);
double RecursiveNestedNewtonForm_getCoef(vector<double>&a, vector<double>&c, double x);
double LagrangeInterpolation(vector<double> &x_i, vector<double>&f_x, double x);
#endif //MN_T9_SOLVERS_H
