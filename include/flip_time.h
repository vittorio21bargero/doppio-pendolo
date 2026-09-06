#ifndef FLIP_TIME_H
#define FLIP_TIME_H
#include "integrators.h"
#include "system.h"
#include <fstream>



void RK4Step (double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq);
double flip_time(double* Y,int neq);
int flip_number(double* Y,int neq);

#endif

