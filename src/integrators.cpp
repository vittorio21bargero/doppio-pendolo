#include "integrators.h"
#include "system.h"
#include <cmath>
#include <string>
#include <vector>

using namespace std;

// the integration methods compared in the project
vector<Method> methods = {
    {"Euler",  EulerStep},
    {"RK2",    RK2Step},
    {"RK4",    RK4Step},
    {"Verlet", position_verlet}
};

// EULER (first order)

void EulerStep(double t, double *Y, void(*RHSFunc)( double*, double*), double dt, int neq) {
    double rhs[4];
    RHSFunc( Y, rhs);

    for (int k = 0; k < neq; k++) {
        Y[k] += dt * rhs[k]; // each component moves along the right-hand side evaluated at the start of the step
    }
}

// RUNGE-KUTTA 2 (midpoint)

void RK2Step (double t, double *Y, void(*RHSFunc)( double*, double*), double h, int neq){

    double Y1[neq], k1[neq], k2[neq];

    RHSFunc(Y,k1);

    int i;

    for (i=0;i<neq;i++){
        Y1[i]= Y[i]+0.5*h*k1[i];}   // half step to the midpoint

    RHSFunc(Y1,k2);                 // slope at the midpoint

    for (i=0;i<neq;i++){
        Y[i]+=h*k2[i];}
}

// RUNGE-KUTTA 4

void RK4Step (double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq){

    double Y1[neq], Y2[neq], Y3[neq], k1[neq], k2[neq], k3[neq], k4[neq];
    int i;

    RHSFunc(Y,k1);
    for (i=0;i<neq;i++){
        Y1[i]= Y[i]+0.5*h*k1[i];}

    RHSFunc(Y1,k2);
    for (i=0;i<neq;i++){
        Y2[i]= Y[i]+0.5*h*k2[i];}

    RHSFunc(Y2,k3);
    for (i=0;i<neq;i++){
        Y3[i]= Y[i]+h*k3[i];}

    RHSFunc(Y3,k4);
    for (i=0;i<neq;i++){
        Y[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);}
}

// POSITION VERLET (drift - kick - drift)
// note: for the double pendulum the accelerations also depend on the velocities,
// so the scheme is not exactly symplectic here

void position_verlet(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq){

    int i;
    double R[neq];

    // positions at t + h/2
    for (i=0; i<neq; i+=2) // i+=2 skips the velocities
        Y[i] += 0.5*h*Y[i+1];

    // accelerations evaluated at the midpoint positions (only R[1], R[3] are used)
    RHSFunc(Y,R);

    // full velocity step, then second half of the position step
    for (i=0; i<neq; i+=2) {
        Y[i+1] += h*R[i+1];
        Y[i]   += 0.5*h*Y[i+1];}
}
