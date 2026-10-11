#ifndef FLIP_TIME_H
#define FLIP_TIME_H

double flip_time(double* Y, int neq);  // first time the second arm flips over (-1 if not within T=10)
int flip_number(double* Y, int neq);   // number of flips of the second arm within T=50

#endif
