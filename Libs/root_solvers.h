#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>
using namespace std;

int Bisection(double (*func)(double), double , double , double , double &zero);
int FalsePos(double (*func)(double), double , double , double , double &zero);
int Secant(double (*func)(double), double , double , double , double &zero);
int Newton(double (*func)(double),double (*dfunc)(double), double , double , double , double &zero);
void Bracket(double (*func)(double),int,int& ,double,double,double [],double []);
