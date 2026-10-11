#ifndef INTEGRATORS_H
#define INTEGRATORS_H
#include <string>
#include <vector>

using namespace std;

// one integration step: advances Y (neq components) by h, using the right-hand side RHSFunc
typedef void (*StepFunc)(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);

struct Method { // pairs the name of a method with the method itself
    string name;
    StepFunc func;
};

void EulerStep(double t, double *Y, void(*RHSFunc)(double*, double*), double dt, int neq);
void RK2Step(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);
void RK4Step(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);
void position_verlet(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);

extern vector<Method> methods; // Euler, RK2, RK4, Verlet

#endif
