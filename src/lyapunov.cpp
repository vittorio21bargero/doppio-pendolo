#include "lyapunov.h"
#include "simulate.h"
#include "integrators.h"
#include <cmath>
#include <fstream>

using namespace std;

double lyapunov_max_benettin(double* Y, int neq)
{
    const double eps = 1e-8;     // initial perturbation
    const double dt = 0.005;     // time between two renormalisations
    const int Nstep = 50;        // RK4 sub-steps per interval --> h = dt/Nstep = 1e-4
    const int transient = 200;   // intervals discarded at the start
    const int Niter = 2000;      // --> total time = 2000*0.005 = 10

    double Y_p[neq]; // perturbed trajectory

    for(int i=0;i<neq;i++) Y_p[i] = Y[i];

    // perturbing one direction is enough: the separation aligns with the most unstable direction
    Y_p[0] += eps;

    double d0 = eps; // separation at t=0
    double d, sum = 0.0;

    ofstream nullfile;

    for(int n=0; n<Niter; n++){

        // advance both trajectories by dt
        integrate(Y , neq, dt, Nstep, RK4Step, dummyCallback, nullfile);
        integrate(Y_p, neq, dt, Nstep, RK4Step, dummyCallback, nullfile);

        // distance in phase space
        d = 0.0;
        for(int i=0;i<neq;i++){
            d += (Y_p[i]-Y[i])*(Y_p[i]-Y[i]);}
        d = sqrt(d);

        // avoid log(0)
        if(d < 1e-14) {
            d = 1e-14;}

        // accumulate only after the transient
        if(n>=transient){
            sum+=log(d/d0);}

        // renormalisation (core of Benettin's algorithm): bring the separation back to d0
        // along the current direction
        for(int i=0;i<neq;i++){
            Y_p[i] = Y[i] + (Y_p[i]-Y[i])*(d0/d);}
    }

    return sum / ((Niter-transient)*dt);
}


double ftle(double* Y, int neq){
    // finite-time Lyapunov exponent from a single separation over T, without renormalisation.
    // Kept for comparison: it saturates once the separation becomes O(1), so it is less reliable.

    const double eps = 1e-8;
    const int N = 10000;   // total steps
    double T=10;
    ofstream nullfile;

    double Y_p[4];

    for(int i=0;i<neq;i++){
        Y_p[i] = Y[i];
    }

    Y_p[0] += eps;

    integrate(Y , neq, T, N, RK4Step, dummyCallback, nullfile);
    integrate(Y_p, neq, T, N, RK4Step, dummyCallback, nullfile);

    double d = 0.0;
    for(int i=0;i<neq;i++)
        d += (Y_p[i]-Y[i])*(Y_p[i]-Y[i]);

    d = sqrt(d);

    return log(d/eps) / T;
}
