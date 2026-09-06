#ifndef INTEGRATORS_H
#define INTEGRATORS_H
#include <string>

using namespace std;

struct Method { // struct mi permette di associare il nome del metodo al metodo vero e proprio
    string name;
    void(*func)(double, double*, void(*)(double*, double*), double, int);
};
void EulerStep(double t, double *Y, void(*RHSFunc)(double*, double*), double dt, int neq);
void RK2Step(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);
void RK4Step(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);
void position_verlet(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);

#endif

