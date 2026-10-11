#ifndef SIMULATE_H
#define SIMULATE_H

#include "integrators.h"
#include <fstream>
#include <string>

using namespace std;

// called at every integration step, e.g. to write the trajectory to file
typedef void (*Callback)(double t, double* Y, ofstream& out, int& k, double& err_cum);

void init_conditions(double* Y, double th1, double w1, double th2, double w2);
void dummyCallback(double t, double* Y, ofstream& out, int& k, double& err_cum);
void writeEnergy(double t, double* Y, ofstream& out, int& k, double& err_cum);
void integrate(double* Y, int neq, double tmax, int N, StepFunc stepFunc,
               Callback callback, ofstream& out);
void integrate_linearized(double* Y, int neq, double tmax, int N, StepFunc stepFunc,
                          Callback callback, ofstream& out);
void compute_grid(const string& filename, double (*quantity_function)(double*, int),
                  int Ngrid = 100);

#endif
