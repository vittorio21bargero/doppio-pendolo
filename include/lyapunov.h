#ifndef LYAPUNOV_H
#define LYAPUNOV_H

double lyapunov_max_benettin(double* Y, int neq); // largest Lyapunov exponent (Benettin et al.)
double ftle(double* Y, int neq);                  // finite-time Lyapunov exponent, single run

#endif
