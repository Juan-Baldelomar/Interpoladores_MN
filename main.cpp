#include <iostream>
#include "solvers.h"
#include "Tools.h"

using  namespace  std;

// funciones a usar
double f2(double x){
    return x*x*(5*x-3)-2*x*x*x*x +4*x - 5;
}

double f3(double x){
    return x*x - 5;
}

double D_f2(double x){
    return 2*x*(5*x-3) + 5*x*x - 8*x*x*x +4;
}

double D_f3(double x){
    return 2*x;
}

void ejercicio1_b(){
    vector<double> x_, f_x;
    x_.assign(4, 0.0);
    f_x.assign(4, 0.0);

    x_[0] = 1981; f_x[0] = 68.3329;
    x_[1] = 1991; f_x[1] = 84.6421;
    x_[2] = 2001; f_x[2] = 102.8737;
    x_[3] = 2011; f_x[3] = 121.0193;

    double x = NDD(x_, f_x, 2006);
    cout << x << endl;
}

void ejercicio2_a(){
    vector<double> x_, f_x, D_f;
    x_.assign(4, 0);
    f_x.assign(4, 0);
    D_f.assign(4, 0);

    x_[0] = -1; f_x[0] = 2; D_f[0] = 2;
    x_[1] = 0; f_x[1] = 2; D_f[1] = 0;
    x_[2] = 1; f_x[2] = 2; D_f[2] = 2;
    x_[3] = 2; f_x[3] = 26; D_f[3] = 68;

    double p = HermiteInterpolation(x_, f_x, D_f, 0.5);
    cout << p << endl;
}

void ejercicio3_a(){
    vector<double> x_, f_x, D_f;
    x_.assign(5, 0);
    f_x.assign(5, 0);

    x_[0] = 1; f_x[0] = 0;
    x_[1] = 1.25; f_x[1] = 0.223144;
    x_[2] = 1.5; f_x[2] = 0.405465;
    x_[3] = 1.75; f_x[3] = 0.559616;
    x_[4] = 2; f_x[4] = 0.693147;

    double p = GregoryNewton(x_, f_x, 1.2);
    cout << p << endl;
}

void ejercicioNDD_sin(){
    int n = 30;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -M_PIf128;
    double b = M_PIf128;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = sin(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = NDD(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = sin(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    //WriteMatrix(points_T, "Output10/sin_NDD_20.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "sin(x)" <<  setw(12) << "NDD" << endl;
    cout << points_T << endl;
}

void ejercicioNDD_f2(){
    int n = 20;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -2;
    double b = 4;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f2(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = NDD(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f2(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f2_NDD_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)" <<  setw(12) << "NDD" << endl;
    cout << points_T << endl;
}

void ejercicioNDD_f3(){
    int n = 20;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -10;
    double b = 10;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f3(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = NDD(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f3(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f3_NDD_20.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)" <<  setw(12) << "NDD" << endl;
    cout << points_T << endl;
}




void ejercicioGN_sin(){
    int n = 10;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -M_PIf128;
    double b = M_PIf128;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = sin(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = GregoryNewton(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = sin(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/sin_GN_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "sin(x)" <<  setw(12) << "GN" << endl;
    cout << points_T << endl;
}

void ejercicioGN_f2(){
    int n = 20;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -2;
    double b = 4;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f2(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = GregoryNewton(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f2(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f2_GN_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)" <<  setw(12) << "GN" << endl;
    cout << points_T << endl;
}

void ejercicioGN_f3(){
    int n = 3;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -10;
    double b = 10;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f3(points[0][i]);
    }

    // NDD calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = GregoryNewton(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f3(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f3_GN_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f3(x)" <<  setw(12) << "GN" << endl;
    cout << points_T << endl;
}


/*  ******************************************************************* HERMITONIAN INTERPOLATION ************************************************************************ */

void ejercicioHer_sin(){
    int n = 10;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(4, vector<double>(n, 0.0));
    double a = -M_PIf128;
    double b = M_PIf128;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = sin(points[0][i]);
        points[2][i] = cos(points[0][i]);
    }

    // Hermite calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[3][i] = HermiteInterpolation(points[0], points[1], points[2], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = sin(x);
        points[2][i] = cos(x);
        //error
        error+= fabs(points[1][i] - points[3][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/sin_Her_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "sin(x)"  << setw(12) << "cos(x)"<<  setw(12) << "Her" << endl;
    cout << points_T << endl;
}

void ejercicioHer_f2(){
    int n = 20;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(4, vector<double>(n, 0.0));
    double a = -2;
    double b = 4;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f2(points[0][i]);
        points[2][i] = D_f2(points[0][i]);
    }

    // Hermite calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[3][i] = HermiteInterpolation(points[0], points[1], points[2], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f2(x);
        points[2][i] = D_f2(x);
        //error
        error+= fabs(points[1][i] - points[3][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f2_Her_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)"  << setw(12) << "D_f2(x)"<<  setw(12) << "Her" << endl;
    cout << points_T << endl;
}

void ejercicioHer_f3(){
    int n = 20;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(4, vector<double>(n, 0.0));
    double a = -10;
    double b = 10;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f3(points[0][i]);
        points[2][i] = D_f3(points[0][i]);
    }

    // Hermite calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[3][i] = HermiteInterpolation(points[0], points[1], points[2], x);
    }

    // update data to see results original function vs lagrange in mid points
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f3(x);
        points[2][i] = D_f3(x);
        //error
        error+= fabs(points[1][i] - points[3][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Output10/f3_Her_10.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f3(x)"  << setw(12) << "D_f3(x)"<<  setw(12) << "Her" << endl;
    cout << points_T << endl;
}

void ej_splines(){
    vector<double> x_, f_;
    x_.assign(4, 0);
    f_.assign(4, 0);

    x_[0] = 0; f_[0] = 1;
    x_[1] = 1; f_[1] = -8;
    x_[2] = 2; f_[2] = -30;
    x_[3] = 3; f_[3] = -59;

    splines s(x_, f_);
    cout << s.P(2.6) << endl;

}


int main(){
    ej_splines();
    return 0;
}