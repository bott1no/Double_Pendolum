#include "ode.h"

using namespace std;


void RK4Step(double t, double *Y, void (*RHS_Func)
(double, double *, double *), int Neq, double h)
{
  double Y1[64], k1[64], k2[64], k3[64], k4[64];
  
  RHS_Func(t,Y,k1); //k1
  
  for (int i=0;i<Neq;i++)
  {
    Y1[i] = Y[i]+0.5*h*k1[i];  //k2
  }
  
  RHS_Func(t+0.5*h,Y1,k2);
  
  for (int i=0;i<Neq;i++)
  {
    Y1[i] = Y[i]+0.5*h*k2[i]; //k3
  }
  
  RHS_Func(t+0.5*h,Y1,k3);
  
  for (int i=0;i<Neq;i++)
  {
    Y1[i] = Y[i]+h*k3[i]; //k4
  }
  
  RHS_Func(t+h,Y1,k4);
  
  for (int i=0;i<Neq;i++)
  {
    Y[i] += h*(k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i])/6.;
  }
}


void RK2Step(double t, double *Y, void (*RHS_Func)(double, double *, double *), int Neq, double h)
{
  double Y1[256], k1[256], k2[256];
  RHS_Func(t,Y,k1);
  
  for (int i=0;i<Neq;i++)
  {
    Y1[i] = Y[i]+0.5*h*k1[i];
  }
  
  RHS_Func(t+0.5*h,Y1,k2);
  
  for (int i=0;i<Neq;i++)
  {
    Y[i] += h*k2[i];
  }
}


void EulerStep (double t, double *Y, void (*RHSFunc)(double, double *, double *), int neq, double dt)
{
int k;
double rhs[256]; //large enough
RHSFunc (t, Y, rhs);
  for (k = 0; k < neq; k++)
  {
    Y[k] += dt*rhs[k];
    //cout<<Y[k]<<endl;
  } // Update solution array
}



