#include "system.h"
#include <cmath>

void dYdt(double *Y, double *R) {
    double m1=1, m2=1, L1=1, L2=1, g=9.81;
    double th1 = Y[0], w1 = Y[1], th2 = Y[2], w2 = Y[3];
    double delta = th1 - th2;
    double den = (2*m1+m2 - m2*cos(2*delta));

    R[0] = w1;
    R[2] = w2;
    R[1] = (-g*(2*m1+m2)*sin(th1) - m2*g*sin(th1-2*th2) - 2*sin(delta)*m2*(w2*w2*L2 + w1*w1*L1*cos(delta))) / (L1*den);
    R[3] = (2*sin(delta)*(w1*w1*L1*(m1+m2) + g*(m1+m2)*cos(th1) + w2*w2*L2*m2*cos(delta))) / (L2*den);
}

void dYdt_linearized(double *Y, double *R) {
    double m1=1, m2=1, L1=1, L2=1, g=9.81;
    double th1 = Y[0], w1 = Y[1], th2 = Y[2], w2 = Y[3];

    R[0] = w1;
    R[2] = w2;
    R[1] = -g*(m1+m2)/(m1*L1) * th1 + (m2*g)/(m1*L1) * th2;
    R[3] =  g*(m1+m2)/(m1*L2) * th1 - g*(m1+m2)/(m1*L2) * th2;
}


// --- exact solution of the linearized system (m1=m2, L1=L2=L), theta1(0)=theta2(0)=theta0, at rest ---
// superposition of the two normal modes w_+ and w_-
double theta1_theoretical(double t, double theta0, double g, double L){
    double w_plus  = sqrt((g/L)*(2.0 + sqrt(2.0)));
    double w_minus = sqrt((g/L)*(2.0 - sqrt(2.0)));
    double c_plus  = 1.0 - 1.0/sqrt(2.0);
    double c_minus = 1.0 + 1.0/sqrt(2.0);
    return 0.5*theta0*( c_plus*cos(w_plus*t) + c_minus*cos(w_minus*t) );
}
