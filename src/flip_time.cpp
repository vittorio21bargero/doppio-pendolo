#include "flip_time.h"
#include "integrators.h"
#include "system.h"
#include <cmath>

using namespace std;

double flip_time(double* Y, int neq) {

    const double dt = 1e-4;      // step chosen from the convergence study
    const double T  = 10.0;      // maximum time
    double t = 0.0;

    while(t < T){
        RK4Step(t, Y, dYdt, dt, neq);   // RK4: the most accurate of the four methods
        t += dt;
        if(fabs(Y[2]) > M_PI){         // flip condition: the second arm goes over the top
            return t;
        }
    }
    return -1.0;   // no flip within T
}


int flip_number(double* Y, int neq) {

    const double dt = 1e-3;
    const double T  = 50.0;
    double t = 0.0;
    int flips = 0;

    // bands of width 2pi centred on 0, 2pi, 4pi, ...: crossing into a new band is a flip
    int k_prev = (int)floor((Y[2] - M_PI) / (2.0*M_PI));

    while (t < T) {
        RK4Step(t, Y, dYdt, dt, neq);
        t += dt;

        int k_curr = (int)floor((Y[2] - M_PI) / (2.0*M_PI));
        if (k_curr != k_prev) {
            flips += 1;
            k_prev = k_curr;
        }
    }
    return flips;
}
