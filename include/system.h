#ifndef SYSTEM_H
#define SYSTEM_H

// State vector: Y = (theta1, omega1, theta2, omega2)
// Parameters: m1 = m2 = 1, L1 = L2 = 1, g = 9.81

void dYdt(double *Y, double *R);              // full nonlinear equations of motion
void dYdt_linearized(double *Y, double *R);   // small-angle approximation
double theta1_theoretical(double t, double theta0, double g, double L);

#endif
