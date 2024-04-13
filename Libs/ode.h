#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

void EulerStep (double , double *, void (*)(double, double *, double *), int , double);
void RK2Step(double , double *  , void (*)(double , double * , double *),int, double);
void RK4Step(double , double *  , void (*)(double , double * , double *),int, double);

