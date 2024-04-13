#include "root_solvers.h"

using namespace std;

int Bisection(double (*func)(double), double a, double b, double tol, double &zero)
{
  int k = 1;
  double fa,fb,xm,fm;
  
  while (fabs(a-b)>tol)
  {
    xm = (a+b)*0.5;
    fm = func(xm);
    fa = func(a);
    fb = func(b);
    
    //cout<<"Bisection():  k = "<< k <<" ;"<<" [a, b] = [" << a;
    //cout<<" , "<< b <<" ];"<<" xm = "<< xm <<" ; dx = "<<fabs(a-b)<<" ; fm = "<< fm<<endl;
    
    if (fa*fm<0) {
      b = xm;
    }else if (fa*fm>0) {
      a = xm;
    } else if (fm == 0) {
      zero = xm;
      return 0;
    } else if (fa == 0) {
      zero = a;
      return 0;
    } else if (fb == 0) {
      zero = b;
      return 0;
    }
  
    k++;
    zero = xm;
  }
  cout<<"Bisection(): x0 = "<< zero <<endl;
  return 0;
}

int FalsePos(double (*func)(double), double a, double b, double tol, double &zero){
  
  int k = 1;
  double xm = 1.;
  double del=1.e3;
  double x0 ;
  double fm,fa,fb;
  
  do{
    
    fa = func(a);
    fb = func(b);
    x0 = xm;
    xm = (a * fb - b * fa) / (fb - fa);
    fm = func(xm);
    
    if (fa*fm < 0.) {    //   P   M    M
      del = fabs(a - xm);
      b   = xm;
    } else if (fa*fm>0.) {   // P  P  M
      del = fabs(xm - b);
      a   = xm;
    } else if (fb == 0.0) {
      zero = b;
      return 0;
    } else if (fa == 0.0) {
      zero = a;
      return 0;
    } else if (fm == 0.0) {
      zero = xm;
      return 0;
    }
    k++;
    zero = xm;
    
  } while (fabs(x0-xm)>tol);
  
  cout<<"FalsePos(): x0 = "<< zero <<endl;
  return 0;
}

int Secant(double (*func)(double), double a, double b, double tol, double &zero){
  
  int k = 1;
  double dx;
  
  double fm,fa,fb;
  fa = func(a);
  fb = func(b);
  
  
  if(fa == 0){
    zero = a;
    return 0;
  } else if(fb == 0)
  {
    zero = b;
  }
  
  while (fabs(b-a)>tol)
  {
    //fa = func(a);
    //fb = func(b);
    dx = fb*(b-a)/(fb-fa);
    a = b;
    //fa = fb;
    b = b - dx;
    fb = func(b);
    k++;
    zero = a;
  }
  
  cout<<"Secant():  k = "<< k <<" ;"<<" xa = " << a <<" ; xb =  "<< b <<" ;"<<" x0 = "<< zero <<endl;
  return 0;
}

int Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, double &zero){
  
  int k = 1;
  double dx,xc;
  double fm,fc;
  xc = (a+b)*0.5;
  //b = mean;

  do{
    
    fc = func(xc);
    dx = fc/dfunc(xc);
    xc -= dx;
    k++;
    
    if (func(xc) == 0.0) {
      break;
    }
    
  } while (fabs(dx)>tol);
  
  cout<<"Newton():  k = "<< k <<" ;"<<" dx = "<< dx <<" ; xc = "<< xc <<endl;
  return 0;
}

void Bracket(double (*func)(double), int N, int& Nr ,double a, double b, double xL[], double xR[]){
  
  double dx = fabs(a - b)/double(N);
  int k = 0;
  
  for (int i = 0; i<N; i++) {
    
    if (func(a)*func(a+dx)<0) {
      xL[k] = a;
      xR[k] = a + dx;
      Nr+=1;
      k++;
    }
    a += dx;
  }
  
}
