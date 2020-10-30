#include <iostream>
#include "solvers.h"
#include "Tools.h"

using  namespace  std;

// funciones a usar
double f2_x(double x){
    return x*x*(5*x-3)-2*x*x*x*x +4*x - 5;
}

double f3_x(double x){
    return x*x - 5;
}

void ejercicio1(){
    int n = 5;
    vector<double> a;
    vector<double> c;

    a.assign(n, 0.0);
    c.assign(n-1, 0.0);

    //init a
    for (int i = 0; i<n; i++)
        a[i]= i;
    // init c
    for (int i = 0; i<n-1; i++)
        c[i] = i;

    double x = 4;
    double r = RecursiveNestedNewtonForm_getCoef(a, c, x);
    cout << "points: " << a << endl;
    cout << "f(x):" << r << endl;
}

void ejercicio2_sin(){
    int n = 100;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -M_PIf128;
    double b = M_PIf128;
    double delta = b-a;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = sin(points[0][i]);
    }

    for (int i = 0; i<n; i++){
        points[2][i] = LagrangeInterpolation(points[0], points[1], points[0][i]);
        error += fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/sin_100.dat");

    // print data in console
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "sin(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;
}

void ejercicio2_sin_diffPoints(){
    int n = 100;
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

    // lagrange calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = LagrangeInterpolation(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = sin(x);
        //error
        error+= fabs(points[1][i] - points[2][i]);
    }

    error /= n;


    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/sin_50_diff.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "sin(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;

}

void ejercicio2_f2(){
    int n = 100;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -2;
    double b = 4;
    double delta = b-a;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f2_x(points[0][i]);
    }

    for (int i = 0; i<n; i++){
        points[2][i] = LagrangeInterpolation(points[0], points[1], points[0][i]);
        error += fabs(points[1][i] - points[2][i]);
    }
    error /= n;


    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/f2_100.dat");

    // print data in console
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;
}

void ejercicio2_f2_diffPoints(){
    int n = 5;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));

    //points to use
    double a = -2;
    double b = 4;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;

    // calculate f(x)
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f2_x(points[0][i]);
    }

    // lagrange calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = LagrangeInterpolation(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f2_x(x);
        // error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    //output
    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/f2_5_diff.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f2(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;

}

void ejercicio2_f3(){
    int n = 100;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));
    double a = -10;
    double b = 10;
    double delta = b-a;
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f3_x(points[0][i]);
    }

    for (int i = 0; i<n; i++){
        points[2][i] = LagrangeInterpolation(points[0], points[1], points[0][i]);
        error += fabs(points[1][i] - points[2][i]);
    }
    error /= n;


    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/f3_100.dat");

    // print data in console
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f3(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;
}

void ejercicio2_f3_diffPoints(){
    int n = 3;
    double error = 0;
    vector<vector<double>> points, points_T;
    points.assign(3, vector<double>(n, 0.0));

    //points to use
    double a = -10;
    double b = 10;
    double delta = b-a;
    double midStep = (delta/(n-1))/2;

    // calculate f(x)
    for (int i = 0; i<n; i++){
        points[0][i] = i*(delta/(n-1)) + a;
        points[1][i] = f3_x(points[0][i]);
    }

    // lagrange calculation
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[2][i] = LagrangeInterpolation(points[0], points[1], x);
    }

    // update data to see results original function vs lagrange
    for (int i = 0; i<n; i++){
        double x = i*(delta/(n-1)) + a+ midStep;
        points[0][i] = x;
        points[1][i] = f3_x(x);
        // error
        error+= fabs(points[1][i] - points[2][i]);
    }
    error /= n;

    //output
    TransposeMatrix(points, points_T);
    WriteMatrix(points_T, "Input/f3_3_diff.dat");
    cout << "Error Promedio: " << error << endl;
    cout << setw(12) <<  "x"  << setw(12) << "f3(x)" <<  setw(12) << "Lagrange" << endl;
    cout << points_T << endl;

}


void pruebaLagrange(){
    vector<double> x_i;
    vector<double> f_x;
    x_i.push_back(1981); x_i.push_back(1991); x_i.push_back(2001); x_i.push_back(2011);
    f_x.push_back(68.3329); f_x.push_back(84.6421); f_x.push_back(102.8737); f_x.push_back(121.0193);
    cout << LagrangeInterpolation(x_i, f_x, 2006);
}


int main() {
    cout << "******************************************************** ALGORITMO RECURSIVO ********************************************************" <<endl;
    ejercicio1();
    cout << "******************************************************** x*x*(5*x-3)-2*x*x*x*x +4*x - 5 ********************************************************" <<endl;
    ejercicio2_f2();
    cout << "******************************************************** x*x*(5*x-3)-2*x*x*x*x +4*x - 5 ********************************************************" <<endl;
    ejercicio2_f2_diffPoints();
    cout << "******************************************************** x^2-10 ********************************************************" <<endl;
    ejercicio2_f3();
    cout << "******************************************************** x^10 - 10 ********************************************************" <<endl;
    ejercicio2_f3_diffPoints();
    cout << "******************************************************** sin (x) ********************************************************" <<endl;
    ejercicio2_sin();
    cout << "******************************************************** sin(x) ********************************************************" <<endl;
    ejercicio2_sin_diffPoints();
    return 0;
}
