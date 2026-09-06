#ifndef SYSTEM_H
#define SYSTEM_H

void dYdt(double *Y, double *R);
void dYdt_linearized(double *Y, double *R);
double theta1_theorical(double t, double theta0, double g, double L);

#endif

