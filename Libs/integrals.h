#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double Simpson (double(*)(double), double, double, int);
double GaussQuadratureRule (double (*)(double), double , double , int , int);
double Trapezoidal (double (*)(double), double , double , int );
double Rectangular (double (*)(double), double , double , int );
