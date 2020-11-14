//
// Created by juan on 23/10/20.
//

#ifndef MN_T9_SOLVERS_H
#define MN_T9_SOLVERS_H

#include "math.h"
#include <vector>
#include "Tools.h"

using namespace std;

class splines{
private:
    vector<double>M;
    vector<double>x_;
    vector<double>f_;
public:
    splines(vector<double> &x_, vector<double> &f_);
    double P(double x);
};


//utilities
void TransposeMatrix(vector<vector<double>>& A, vector<vector<double>>&B);

double RecursiveNestedNewtonForm(vector<double>&a, vector<double>&c, double x);
double RecursiveNestedNewtonForm_getCoef(vector<double>&a, vector<double>&c, double x);
double LagrangeInterpolation(vector<double> &x_i, vector<double>&f_x, double x);

double NDD(vector<double>&x_, vector<double>&f_x, double x);
double HermiteInterpolation(vector<double> &x_, vector<double>&f_x, vector<double>&Dx_f, double x);
double GregoryNewton(vector<double>&x_, vector<double>& f_x, double x);
#endif //MN_T9_SOLVERS_H
