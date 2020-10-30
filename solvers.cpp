//
// Created by juan on 23/10/20.
//

#include "solvers.h"

void TransposeMatrix(vector<vector<double>>& A, vector<vector<double>>&B){
    int n = A.size();
    int m = A[0].size();

    B.assign(m, vector<double>(n, 0.0));
    for (int i = 0; i<n; i++)
        for (int j= 0; j<m; j++)
            B[j][i] = A[i][j];
}

double l_i(vector<double>&x_i, double x, int i){
    int n = x_i.size();
    double acc = 1;
    for (int k = 0; k<n; k++){
        if (k==i)
            continue;
        acc*= (x-x_i[k])/(x_i[i] - x_i[k]);
    }
    return acc;
}

// does not change coeficients a
double RecursiveNestedNewtonForm(vector<double>&a, vector<double>&c, double x){
    int n = a.size();
    double inner_exp = a[n-1];                          // a'_n = a_n
    for (int i = n-1; i>1; i--){
        inner_exp = a[i-1]+ inner_exp*(x-c[i-1]);       // a'_{n-1} = a_{n-1} + (x-c_n)a'_n
    }
    return a[0] + inner_exp*(x-c[0]);
}

double RecursiveNestedNewtonForm_getCoef(vector<double>&a, vector<double>&c, double x){
    int n = a.size();
    // inner_exp stored in a[i]
    for (int i = n-1; i>1; i--){
        a[i-1] = a[i-1]+ a[i]*(x-c[i-1]);       // a'_{i-1} = a_{i-1} + (x-c_n)a'_i
    }
    a[0] =a[0] + a[1]*(x-c[0]);
    return a[0];
}

double LagrangeInterpolation(vector<double> &x_i, vector<double>&f_x, double x){
    int n = x_i.size();
    double p_x = 0;
    for (int i = 0; i<n; i++){
        p_x += f_x[i]*l_i(x_i, x, i);
    }
    return p_x;
}