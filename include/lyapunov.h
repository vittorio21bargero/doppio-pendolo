#ifndef LYAPUNOV_H
#define LYAPUNOV_H

#include <fstream>
 using namespace std;






void integrate(double* Y, int neq, double tmax, int N,
               void(*stepFunc)(double, double*, void(*)(double*, double*), double, int),
               void(*callback)(double, double*,ofstream&,int&,double&),
               ofstream& out);
double lyapunov_max_benettin(double* Y, int neq);
double ftle(double* Y, int neq);

#endif
