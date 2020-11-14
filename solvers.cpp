//
// Created by juan on 23/10/20.
//

#include "solvers.h"
#include "iostream"

using namespace std;

void TransposeMatrix(vector<vector<double>>& A, vector<vector<double>>&B){
    int n = A.size();
    int m = A[0].size();

    B.assign(m, vector<double>(n, 0.0));
    for (int i = 0; i<n; i++)
        for (int j= 0; j<m; j++)
            B[j][i] = A[i][j];
}

void Triangular_Superior(vector<vector<double> > &U, vector<double> &b, vector<double> &x) {
    int n = U.size();
    double X[n];

    for (int i = n - 1; i >= 0; i--) {
        double acc = 0;
        for (int j = n - 1; j > i; j--) {
            acc += U[i][j] * x[j];
        }
        if (U[i][i] == 0)
            cout << "Triangular Sup  Msg WARNING: A_[" << i << ", " << i << "] es nulo, se dividira por 0" << endl;

        x[i] = (b[i] - acc) / U[i][i];
    }
}

void Partial_pivot(vector<vector<double> > &A, int k, vector<double> &b) {
    int n = A.size();
    int i_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k + 1; i < n; i++) {
        if (fabs(A[i][k]) > max) {
            i_max = i;
            max = fabs(A[i][k]);
        }
    }

    double temp;
    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }
}

void Eliminacion_Gaussiana(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {

    int n = A.size();
    //indice k indica la fila con la cual estamos transformando a una Triangular Superior
    for (int k = 0; k < n - 1; k++) {

        if (A[k][k] == 0)
            Partial_pivot(A, k, b);

        if (A[k][k] == 0)
            cout << "GEM  Msg WARNING: A_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }
    }
    Triangular_Superior(A, b, x);
}

int factorial(int n){
    int acc = 1;
    for (int i = 1; i<= n; i++)
        acc*= i;

    return acc;
}

double exp(double x, int n){
    double acc = 1;
    for (int i = 0; i<n; i++)
        acc*= x;
    return acc;
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

double Dx_l_i(vector<double>&x_, double x, int i){
    int n = x_.size();
    double num = 0;
    double den = 1;
    for (int k = 0; k<n; k++){
        if (k==i)
            continue;
        den*= (x_[i] - x_[k]);
    }

    for (int k = 0; k<n; k++){
        if (k==i)
            continue;
        double mult = 1;
        for (int j = 0; j<n; j++){
            if (j==k || j==i)
                continue;
            mult*= (x-x_[j]);
        }
        num += mult;
    }
    return num/den;
}

// x_ is the dataset from which we want to get the forward Difference Table
// FD_Table must have not been assigned a size, must be empty
void Forward_Difference(vector<double>&x_, vector<vector<double>>&FD_Table){
    int n = x_.size();

    FD_Table.push_back(vector<double>(n, 0.0));
    for(int i = 0; i<n; i++){
        FD_Table[0][i] = x_[i];
    }

    for (int k = 1; k<n; k++){
        FD_Table.push_back(vector<double>(n-k, 0.0));
        int size = FD_Table[k].size();
        for (int i = 0; i<size; i++){
            FD_Table[k][i] = (FD_Table[k-1][i+1] - FD_Table[k-1][i]);
        }
    }
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

double x_prod(vector<double> &x_i, double x, int k){
    if (k==0)
        return 1;

    double acc = 1;
    for (int i = 0; i<k; i++){
        acc*= x-x_i[i];
    }
    return acc;
}

/* ****************************************************************************** PRACTICA ACTUAL ****************************************************************/

double NDD(vector<double>&x_, vector<double>&f_x, double x){
    int n = x_.size();
    vector<vector<double>> DDTable;
    vector<double> a;

    DDTable.push_back(vector<double>(n, 0.0));

    // init first column
    for (int i = 0; i<n; i++){
        DDTable[0][i] = f_x[i];
    }
    a.push_back(DDTable[0][0]);                                 //a_0

    // calculate the other columns
    // Notice this is a triangular matrix
    for (int k = n-1, delta_X = 1; k>0; k--, delta_X++){
        DDTable.push_back(vector<double>(k, 0.0));
        int dim = n-(k+1);                                      //dimension to calculate new dimension
        int size = DDTable[dim].size();
        for (int i = 0; i<size-1; i++){
            DDTable[dim+1][i] = (DDTable[dim][i+1] - DDTable[dim][i])/(x_[i+delta_X] - x_[i]);
        }
        a.push_back(DDTable[dim+1][0]);                         //a_i
    }
    return RecursiveNestedNewtonForm_getCoef(a, x_, x);
}

double HermiteInterpolation(vector<double> &x_, vector<double>&f_x, vector<double>&Dx_f, double x) {
    int n = x_.size();
    double acc = 0;
    for (int i = 0; i<n; i++){
        double l_ = l_i(x_, x, i);
        double D_l_ = Dx_l_i(x_, x_[i], i);
        double u_i = (-2 * D_l_ * x + 1 + 2 * x_[i] * D_l_) * l_ * l_ * f_x[i];
        double v_i = (x-x_[i]) * l_ * l_ * Dx_f[i];
        acc += u_i + v_i;
    }
    return acc;
}

double GregoryNewton(vector<double>&x_, vector<double>& f_x, double x){
    int n = x_.size();
    double h = x_[1] - x_[0];
    vector<vector<double>> FD_Table;
    vector<double>a;
    Forward_Difference(f_x, FD_Table);

    a.push_back(f_x[0]);
    for (int k = 1; k<n; k++){
        double quotient = 1/(factorial(k)*exp(h, k));
        double a_i= quotient * FD_Table[k][0];
        a.push_back(a_i);
    }

    return RecursiveNestedNewtonForm(a, x_, x);
}


// **************************************************************************** SPLINES CLASS ****************************************************************************


splines::splines(vector<double> &x_, vector<double> &f_) {
    int n = x_.size();

    M.assign(n-2, 0);
    this->x_.assign(n, 0);
    this->f_.assign(n, 0);

    double h = x_[1] - x_[0];

    //build system
    vector<vector<double>> A(n - 2, vector<double>(n - 2, 0));
    vector<double> b(n-2, 0);

    for (int i = 0; i<n; i++){
        this->x_[i] = x_[i];
        this->f_[i] = f_[i];
    }

    for (int i = 0; i<n-2; i++){
        if (i == 0){
            A[i][0] = 4;
            A[i][1] = 1;
        }else if (i == n-3){
            A[i][n - 4] = 1;
            A[i][n - 3] = 4;
        }else{
            A[i][i - 1] = 1;
            A[i][i] = 4;
            A[i][i + 1] = 1;
        }
                                                                        // 6/h^2 * (f_{i+1} -2 f_i + f_{i-1})
        b[i] = (6/(h*h)) * (f_[i+2] - 2*f_[i+1] + f_[i] );              // i -> i+1 porque i inicia desde 0 en el ciclo y deberia iniciar en 1
    }

    // calculate M
    Eliminacion_Gaussiana(A, b, M);
    M.insert(M.begin(), 0);             //M_0 = 0
    M.insert(M.end(), 0);               //M_n = 0
}

double splines::P(double x) {
    int n = x_.size();
    double h = x_[1] - x_[0];
    int i;
    for ( i = 1; i<n-1; i++){
        if (x_[i-1]<= x && x<= x_[i])
            break;
    }

    double a, b, c;

    a = 1/(6*h) * (pow(x_[i] - x, 3) * M[i-1] + pow(x - x_[i-1], 3)* M[i]);
    b = ((x_[i] - x)/h ) * (f_[i-1]- (h*h/6) * M[i-1]);
    c = ((x-x_[i-1])/h ) * (f_[i] - (h*h)/6 * M[i]);

    return a + b + c;
}