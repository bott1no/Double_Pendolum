#include "ode.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

void dYdt (double , double *, double *);

#define STAGE 1

//STAGE 1 == convergency test
//STAGE 2 == gif with 2 systems
//STAGE 3 == finding the configurations that lead to chaos
//STAGE 4 == energy conservation
//STAGE 5 == finding the configurations that lead to flip

int main(){
  
  double t = 0.;
  double tf = 20.;
  double dt = (tf - t)/double(2048);
  double L1 = 1., L2 = 1., g = 9.81;
  int neq = 4;
  double R[256];
  double gg = 9.81;
  
  //int N = 10;
  
#if STAGE == 1 //convergency test
  
  int N = 1024;

  double Y4[256];
  double Y2[256];
 //setting the error at an high value
  double err2 ;//setting the error at an high value
  double err4 ;//setting the error at an high value
  double err1 ;//setting the error at an high value
  L2 = 0.5;
  L1 = 1.;
  double L = L2/double(L1);
  
  ofstream convergency;
  convergency.open("convergency.dat");
  Y4[0] = 0.; //setting BC
  Y4[1] = 0.;
  Y4[2] = 1./100.; //setting a small angle oscillation
  Y4[3] = 1./100.;
  
  Y2[0] = 0.; //setting BC
  Y2[1] = 0.;
  Y2[2] = 1./100.; //setting a small angle oscillation
  Y2[3] = 1./100.;
  double BC1 = Y4[2];
  double BC2 = Y4[3];
  double BC3 = Y4[0];
  double BC4 = Y4[1];
  tf = 10.;
 
  while (N < 400000) {
    t = 0.;
    dt = (tf - 0.)/double(N);
    Y4[0] = 0.; //setting BC
    Y4[1] = 0.;
    Y4[2] = 1./100.; //setting a small angle oscillation
    Y4[3] = 1./100.;
    
    Y2[0] = 0.; //setting BC
    Y2[1] = 0.;
    Y2[2] = 1./100.; //setting a small angle oscillation
    Y2[3] = 1./100.;
    
    for(int l = 0; l<N;l+=1)
    {
      RK4Step(t, Y4, dYdt, neq, dt);
      RK2Step(t, Y2, dYdt, neq, dt);
      t += dt;
    }
    double omegaQuadro = gg/double(L1);
    double exact2 = (BC2 - BC1/(1.-double(L))) * cos(tf/sqrt(double(L))) + BC1 * cos(tf)/(1.-double(L));
   // double exact2 = ((-BC3+omegaQuadro*(BC1+L*BC3))*cos(t/sqrt(double(L)))+sqrt(L)*(-BC4+omegaQuadro*(BC2+L*BC4))*sin(t/sqrt(double(L)))-sqrt(omegaQuadro)*(sqrt(omegaQuadro)*BC1*cos(t*sqrt(omegaQuadro)+BC2*sin(t*sqrt(omegaQuadro)))))/(-1.+L*omegaQuadro);
    
    err4 = fabs(Y4[3] - exact2);
    err2 = fabs(Y2[3] - exact2);
    err1 = fabs(Y2[2] - BC1*cos(tf));
    
    
    convergency << dt << " " << 1 << " " << err2 << " " << err4 << " " << endl;
    cout<<N<<"   "<<Y4[3]<<"   "<<exact2<<"  "<<err4<<endl;
    N *= 2;
  }
  convergency.close();
  
#endif
  
  //--------------------------------
  
#if STAGE == 2 // (grafical) difference between two systems same 
               // initial conditions
               // (file doubledouble.gp)
  double Y2[256], Y4[256];
  
  Y2[0] = 0.; //velocities
  Y2[1] = 0.;
  Y2[2] = 1.; //theta1
  Y2[3] = 3.; //theta2
  
  Y4[0] = Y2[0]; //here we set the same BC for RK2 and RK4
  Y4[1] = Y2[1];
  Y4[2] = Y2[2];
  Y4[3] = Y2[3]+1.e-5;
  
  ofstream rk4;
  ofstream rk2;
  rk4.open("DoublePendolum1.dat"); //the blue pendolum in the gif
  rk2.open("DoublePendolum2.dat"); //the red pendolum in the gif
  
  double sin2 = sin(Y4[2]);
  double cos2 = cos(Y4[2]);
  double sin3 = sin(Y4[3]);
  double cos3 = cos(Y4[3]);
  
  double sin2_2 = sin(Y2[2]);
  double cos2_2 = cos(Y2[2]);
  double sin2_3 = sin(Y2[3]);
  double cos2_3 = cos(Y2[3]);
  
  double x1_4 = L1 * sin2;  //definition of the cartesian coordinates
  double y1_4 = - L1 * cos2;
  double x2_4 = x1_4 + L2 * sin3;
  double y2_4 = y1_4 - L2 * cos3;
  
  double x1_2 = L1 * sin2_2;
  double y1_2 = - L1 * cos2_2;
  double x2_2 = x1_2 + L2 * sin2_3;
  double y2_2 = y1_2 - L2 * cos2_3;
  
  t = 0.;

  while (t<=tf)
  {
    sin2 = sin(Y4[2]);
    cos2 = cos(Y4[2]);
    sin3 = sin(Y4[3]);
    cos3 = cos(Y4[3]);
    
    sin2_2 = sin(Y2[2]);
    cos2_2 = cos(Y2[2]);
    sin2_3 = sin(Y2[3]);
    cos2_3 = cos(Y2[3]);
    
    x1_4 = L1 * sin2;
    y1_4 = - L1 * cos2;
    x2_4 = x1_4 + L2 * sin3;
    y2_4 = y1_4 - L2 * cos3;
    
    x1_2 = L1 * sin2_2;
    y1_2 = - L1 * cos2_2;
    x2_2 = x1_2 + L2 * sin2_3;
    y2_2 = y1_2 - L2 * cos2_3;
    
    rk4 << 0 << " " << Y4[0] << " " << Y4[1] <<" "<< x1_4 << " " << y1_4 << " " << x2_4 << " " << y2_4 << " " << endl;
    
    rk2 << 0 << " " << Y2[0] << " " << Y2[1] <<" "<< x1_2 << " " << y1_2 << " " << x2_2 << " " << y2_2 << " " << endl;
  
    RK4Step(t, Y4, dYdt, neq, dt);
    RK2Step(t, Y2, dYdt, neq, dt);
    t += dt;
  }
  rk2.close();
  rk4.close();
#endif
  
#if STAGE == 3  //(careful, long run)
  //showing chaos with slightly different initial conditions
  
  double Y2[256]; //perturbated
  double Y1[256]; //non perturbated
  double perturbation = 1.e-5;
  double sumSquareDiff = 0.;
  double t0 = 0.;
  tf = 20.;
  dt = (tf - t0)/double(2.e3);
  double LyapunovNorm;

  ofstream chaos;
  chaos.open("chaos.dat");
  
  for (double n = 0.01; n<2*M_PI; n+=0.01) {
    for (double m = 0.01; m<2*M_PI; m+=0.01) {
      
      Y1[0] = 0.;
      Y1[1] = 0.;
      Y1[2] = M_PI*2.-n; //theta1
      Y1[3] = M_PI*2.-m; //theta2
      double InY1 = Y1[2];
      double InY2 = Y1[3];
      
      Y2[0] = Y1[0];
      Y2[1] = Y1[1];
      Y2[2] = Y1[2];
      Y2[3] = Y1[3] + perturbation;
    
      t = 0.;
      double sumSquareDiff = 0.;
      while (t <= tf)
      {
        RK4Step(t, Y1, dYdt, neq, dt);
        RK4Step(t, Y2, dYdt, neq, dt);
        t += dt;
        
        for (int i=0; i<4; i++) {
          double diff = fabs(Y1[i] - Y2[i]);
          sumSquareDiff += diff*diff;
        }
        double Lyapunov = log(sqrt(sumSquareDiff));
        LyapunovNorm = Lyapunov/(t-t0);
      }
      if (LyapunovNorm > 0.) {
        chaos << InY1 << " " << InY2 << " " << endl;
        
        /*cout<<"Lyapunov: " << LyapunovNorm <<endl;
        cout<<"starting angles: Pi/" << m <<endl;
        cout<<"Lyapunov exponent: " << LyapunovNorm <<endl<<endl<<endl;*/
     }
    }
    
  }
  chaos.close();
#endif
  
  //--------------------------------
  
#if STAGE == 4 // just the one pendulum --> animated gif 
               // and energy considerations
  
  double Y4[256];
  
  Y4[0] = 0.; //setting BC
  Y4[1] = 0.;
  Y4[2] = 4.;
  Y4[3] = 1.;
  L1=1.;
  L2=1.;
  
  ofstream rk4;
  rk4.open("DoublePendolum1.dat");
  
  ofstream energy;
  energy.open("energy.dat");
  
  double sin2 = sin(Y4[2]);
  double cos2 = cos(Y4[2]);
  double sin3 = sin(Y4[3]);
  double cos3 = cos(Y4[3]);
  
  double x1_4 = L1 * sin2;  //definition of the cartesian coordinates
  double y1_4 = - L1 * cos2;
  double x2_4 = x1_4 + L2 * sin3;
  double y2_4 = y1_4 - L2 * cos3;
  
  double vx1_4 = L1 * cos2 * Y4[0];
  double vy1_4 = L1 * sin2 * Y4[0];
  double vx2_4 = vx1_4 + L2 * cos3 * Y4[1];
  double vy2_4 = vy1_4 + L2 * sin3 * Y4[1];
  
  double Ereal = gg * (y1_4 + y2_4) + 0.5 * (vx1_4 * vx1_4 + vx2_4 * vx2_4 + vy1_4 * vy1_4 + vy2_4 * vy2_4);
                 //evaluate the exact energy with the initial
                 //conditions (potential energy)
  double E4;
  
  t = 0.;
  int N=20000;
  for (int i=0; i<N; i++) {
    sin2 = sin(Y4[2]);
    cos2 = cos(Y4[2]);
    sin3 = sin(Y4[3]);
    cos3 = cos(Y4[3]);
    
    x1_4 = L1 * sin2;
    y1_4 = - L1 * cos2;
    x2_4 = x1_4 + L2 * sin3;
    y2_4 = y1_4 - L2 * cos3;
    
    vx1_4 = L1 * cos2 * Y4[0];
    vy1_4 = L1 * sin2 * Y4[0];
    vx2_4 = vx1_4 + L2 * cos3 * Y4[1];
    vy2_4 = vy1_4 + L2 * sin3 * Y4[1];
    
    E4 = gg * (y1_4 + y2_4) + 0.5 * (vx1_4 * vx1_4 + vx2_4 * vx2_4 + vy1_4 * vy1_4 + vy2_4 * vy2_4);
    
    rk4 << 0 << " " << Y4[0] << " " << Y4[1] << " "<< x1_4 << " " << y1_4 << " " << x2_4 << " " << y2_4  << " " << endl;
    
    energy << t << " " << 1 << " " << E4 << " " << Ereal << " " << fabs((Ereal - E4)/Ereal) << " " << 1 << " " << endl;
    
    RK4Step(t, Y4, dYdt, neq, dt);
    t += dt;
  }
  cout<<Y4[3]<<endl;
  rk4.close();
  energy.close();
#endif
  
#if STAGE == 5 //checking the flip

  double Y1[256];
  double t0 = 0.;
  tf = 20.;
  double N=2.e3;
  dt = (tf - t0)/double(N);
  int count=0;
  ofstream flip;
  flip.open("flip.dat");
  
  
  for (double n = 0.01; n<2*M_PI; n+=0.01) {
    for (double m = 0.01; m<2*M_PI; m+=0.01) {
      
      Y1[0] = 0.;
      Y1[1] = 0.;
      Y1[2] = M_PI*2.-n; //theta1
      Y1[3] = M_PI*2.-m; //theta2
      double InY1 = Y1[2];
      double InY2 = Y1[3];
      t = 0.;
      count = 0;
      for (int k=0; k<N; k++) {
        double pre_theta1 = Y1[2];
        double pre_theta2 = Y1[3];
      
        RK4Step(t, Y1, dYdt, neq, dt);
        t += dt;
        
        if ((pre_theta1 <= M_PI && Y1[2] > M_PI) || (pre_theta1 >= M_PI && Y1[2] < M_PI)) {
            count++;
          
        }else if ((pre_theta2 <= M_PI && Y1[3] > M_PI) || (pre_theta2 >= M_PI && Y1[3] < M_PI)) {
          count++;
          
        }
      }
      
     
        flip << InY1 << " " << InY2 << " " << count << " " << endl;
      if (count>=10) {
        cout<<count<<endl;
      }
     
    }
  }
  flip.close();
  
#endif
#if STAGE == 6 // just the one pendulum --> animated gif
               // and energy considerations

  int N = 8;

  double Y4[256];
 //setting the error at an high value
  double err2 ;//setting the error at an high value
  L2 = 0.5;
  L1 = 1.;
  double L = L2/double(L1);
  
  ofstream convergency;
  convergency.open("convergency.dat");
  Y4[0] = 0.; //setting BC
  Y4[1] = 0.;
  Y4[2] = 3; //setting a small angle oscillation
  Y4[3] = 2.;
  
  tf = 10.;
  
  double E4;
  
  t = 0.;
  
  double sin2 = sin(Y4[2]);
  double cos2 = cos(Y4[2]);
  double sin3 = sin(Y4[3]);
  double cos3 = cos(Y4[3]);
  
  double x1_4 = L1 * sin2;  //definition of the cartesian coordinates
  double y1_4 = - L1 * cos2;
  double x2_4 = x1_4 + L2 * sin3;
  double y2_4 = y1_4 - L2 * cos3;
  
  double vx1_4 = L1 * cos2 * Y4[0];
  double vy1_4 = L1 * sin2 * Y4[0];
  double vx2_4 = vx1_4 + L2 * cos3 * Y4[1];
  double vy2_4 = vy1_4 + L2 * sin3 * Y4[1];
  double Ereal = gg * (y1_4 + y2_4) + 0.5 * (vx1_4 * vx1_4 + vx2_4 * vx2_4 + vy1_4 * vy1_4 + vy2_4 * vy2_4);
  
  while (N < 40000) {
    t = 0.;
    dt = (tf - 0.)/double(N);
    Y4[0] = 0.; //setting BC
    Y4[1] = 0.;
    Y4[2] = 3; //setting a small angle oscillation
    Y4[3] = 2;
    sin2 = sin(Y4[2]);
    cos2 = cos(Y4[2]);
    sin3 = sin(Y4[3]);
    cos3 = cos(Y4[3]);
    
    x1_4 = L1 * sin2;  //definition of the cartesian coordinates
    y1_4 = - L1 * cos2;
    x2_4 = x1_4 + L2 * sin3;
    y2_4 = y1_4 - L2 * cos3;
    
    vx1_4 = L1 * cos2 * Y4[0];
    vy1_4 = L1 * sin2 * Y4[0];
    vx2_4 = vx1_4 + L2 * cos3 * Y4[1];
    vy2_4 = vy1_4 + L2 * sin3 * Y4[1];
      
    for(int l = 0; l<N;l+=1)
     {
       RK4Step(t, Y4, dYdt, neq, dt);
       t += dt;
       sin2 = sin(Y4[2]);
       cos2 = cos(Y4[2]);
       sin3 = sin(Y4[3]);
       cos3 = cos(Y4[3]);
       
       x1_4 = L1 * sin2;  //definition of the cartesian coordinates
       y1_4 = - L1 * cos2;
       x2_4 = x1_4 + L2 * sin3;
       y2_4 = y1_4 - L2 * cos3;
       
       vx1_4 = L1 * cos2 * Y4[0];
       vy1_4 = L1 * sin2 * Y4[0];
       vx2_4 = vx1_4 + L2 * cos3 * Y4[1];
       vy2_4 = vy1_4 + L2 * sin3 * Y4[1];
     }
    
   E4 = gg * (y1_4 + y2_4) + 0.5 * (vx1_4 * vx1_4 + vx2_4 * vx2_4 + vy1_4 * vy1_4 + vy2_4 * vy2_4);
   err2 = fabs(Ereal - E4);
   convergency << dt << " " << 1 << " " << err2 << " " << endl;
   cout<<N<<"   "<<err2<<"   "<<Ereal<<"  "<<E4<<endl;
   N *= 2;
  }
  
  convergency.close();
 
#endif
  
  return 0;
}


void dYdt (double t, double *Y, double *R) //equation of motion
{
  double th1    = Y[2];
  double th2    = Y[3];
  double th1dot = Y[0];
  double th2dot = Y[1];
  double m1 = 1.e6, m2 = 1., L1 = 1., L2 = 0.5, g = 9.81;
  R[2] = Y[0];
  R[3] = Y[1];
  R[0] = (-g * (2. * m1 + m2) * sin(th1) - m2 * g * 
          sin(th1 - 2. * th2) - 2. * m2 * sin(th1 - th2) * (th2dot * th2dot * L2 + th1dot * th1dot * L1 * cos(th1 - th2)))/(L1*(2*m1 + m2 - m2 * cos(2.*th1 - 2.*th2)));
          
  R[1] = (2. * sin(th1 - th2) * (th1dot * th1dot * L1 * (m1 + m2) + g * (m1 + m2) * cos(th1) + th2dot * th2dot * L2 * m2 * cos(th1 - th2)))/(L2*(2*m1 + m2 - m2 * cos(2.*th1 - 2.*th2)));
}

