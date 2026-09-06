#include "integrators.h"
#include "system.h"
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

double flip_time(double* Y, int neq) {

    const double dt = 1e-4;      // passo ottimale
    const double T  = 10.0;      // tmax
    double t = 0.0;

    while(t < T){
        RK4Step(t, Y, dYdt, dt, neq);   // RK4 (metodo migliore)
        t += dt;
        if(fabs(Y[2]) > M_PI){     // condizione di flip
            return t;
        }
    }
    return -1.0;   // nessun flip entro T
}


int flip_number(double* Y, int neq) {

    const double dt = 1e-3;
    const double T  = 50.0;
    double t = 0.0;
    int flips = 0;

    // devo creare fascie di 2pi (traslando di pi) per rimanere in mod 2pi
    int k_prev = (int)floor((Y[2] - M_PI) / (2.0*M_PI));

    while (t < T) {
        RK4Step(t, Y, dYdt, dt, neq);
        t += dt;

        // se mi trovo in una fascia diversa da prima allora il pendolo ha flippato
        int k_curr = (int)floor((Y[2] - M_PI) / (2.0*M_PI));
        if (k_curr != k_prev) {
            flips += 1;
            k_prev = k_curr;
        }
    }
    return flips;
}
