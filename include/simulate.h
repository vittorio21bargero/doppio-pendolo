#ifndef SIMULATE_H
#define SIMULATE_H

#include <fstream>
#include <string>


using namespace std;


void init_conditions(double* Y, double th1, double w1, double th2, double w2);
void dummyCallback(double t, double* Y, ofstream& out, int& k,double& err_cum);
void writeData(double t, double* Y,ofstream& out);
void writeData1(double t, double* Y,ofstream& out, int& k, double& err_cum);
void integrate(double* Y, int neq, double tmax, int N,
               void(*stepFunc)(double, double*, void(*)(double*, double*), double, int),
               void(*callback)(double, double*,ofstream&,int&,double&),
               ofstream& out);
void integrate_linearized(double* Y, int neq, double tmax, int N,
               void(*stepFunc)(double, double*, void(*)(double*, double*), double, int),
               void(*callback)(double, double*,ofstream&,int&,double&),
               ofstream& out);
void compute_grid(const string& filename,
                  double (*quantity_function)(double*,int));

#endif
