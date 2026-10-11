#ifndef ENERGY_H
#define ENERGY_H

extern double E_0; // energy of the initial condition, reference for the relative error

double E(double* Y);       // total mechanical energy
double E_err(double* Y);   // relative energy error |E - E_0| / E_0

#endif
