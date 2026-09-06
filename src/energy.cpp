#include "energy.h"
#include "integrators.h"
#include "simulate.h"
#include <cmath>

double E_0=25;

double E(double* Y) {
    // energia meccanica del doppio pendolo

    double L1=1, L2=1, m1=1, m2=1, g=9.81;
    double V = (m1+m2)*g*L1*(1 - cos(Y[0])) + m2*g*L2*(1 - cos(Y[2]));
    double T = 0.5*(m1+m2)*L1*L1*Y[1]*Y[1] + 0.5*m2*L2*L2*Y[3]*Y[3] + m2*L1*L2*cos(Y[0]-Y[2])*Y[1]*Y[3];
    return T + V;
}

double E_err(double* Y){
    // errore relativo sull'energia

    return abs((E(Y)-E_0)/E_0);


}
