#include "integrals.h"
using namespace std;

double Simpson (double (*F)(double), double a, double b, int N)
{
  double h = (b-a)/double(N);
  
  double sum = h * (F(a) + F(b));
  
  double x = a + h;
  
  double w = 4.0;
  
  for (int i = 0; i < N-1; i++) {
    
    
    sum = sum + w * F(x) * h;
    
    x = x + h;
    
    w = 6. - w;
    
  }
  return sum/3.;
}


double GaussQuadratureRule (double (*F)(double), double a, double b, int N, int Ng)
{
  double wg[8];
  double xg[8];
  
  if (Ng == 2) {
    
    
    wg[0] = 1.;
    wg[1] = 1.;
      
  
    xg[0] = -1./sqrt(3.);
    xg[1] =  1./sqrt(3.);
  }
  
  
  else if (Ng == 3) {
    
    
    wg[0] = 5./9.;
    wg[1] = 8./9.;
    wg[2] = 5./9.;
    
    xg[0] = -sqrt(3./5.);
    xg[1] = 0.;
    xg[2] =  sqrt(3./5.);
  
  }
  
  else if (Ng == 4) {
    
    
    wg[0] = (18.-sqrt(30.))/36.;
    wg[1] = (18.+sqrt(30.))/36.;
    wg[2] = (18.+sqrt(30.))/36.;
    wg[3] = (18.-sqrt(30.))/36.;
    
    
    xg[2] = sqrt(3./7.-(2./7.)*sqrt(6./5.));
    xg[1] = -sqrt(3./7.-(2./7.)*sqrt(6./5.));
    xg[3] =  sqrt(3./7.+(2./7.)*sqrt(6./5.));;
    xg[0] = -sqrt(3./7.+(2./7.)*sqrt(6./5.));;
    
  }
  
  double h = (b-a)/double(N);
  
  double sum = 0.0;
  
  for (int i = 0; i < N; i++) { // Loop over intervals
    
    
    double x0 = a + i*h;
    double  x1 = x0 + h; // Define left and right interval endpoints x0 & x1
    
    double sumk = 0.0; // Initialize sum for this interval
    
    double center = (x0 + x1)/2.;
    
    double halfsum = (x1 - x0)/2.;
    
    for (int k = 0; k < Ng; k++)
      
    {
      
      sumk += halfsum * wg[k] * F( halfsum * xg[k] + center) ;  // Apply Gaussian rule to sub-interval
      
       // Add partial sum to total integral
    }
    sum += sumk;
  }
    return sum;
  
  
}

double Trapezoidal (double (*F)(double), double a, double b, int N)
{
  double h = (b-a)/double(N);
  
  double sum = h / 2 * (F(a) + F(b)) ;
  
  double x = a + h;
  
  for (int i = 0; i < N-1; i++) {
    
    sum = sum + F(x) * h;
    
    x = x + h;
  }
return sum;
  
  
}

double Rectangular (double (*F)(double), double a, double b, int N)
{
  double h = (b-a)/double(N);
  
  double sum = 0.;
  
  
  double x = a;
  
  for (int i = 0; i < N; i++) {
    
    
    sum = sum + F(x) * h;
    
    x = x + h;
  }
return sum;
  
  
}
