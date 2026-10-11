#include "simulate.h"
#include "energy.h"
#include "system.h"
#include "integrators.h"
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;


void init_conditions(double* Y, double th1, double w1, double th2, double w2) {
    Y[0] = th1;
    Y[1] = w1;
    Y[2] = th2;
    Y[3] = w2;
}


void dummyCallback(double t, double* Y, ofstream& out, int& k, double& err_cum) {
    // does nothing: used when only the final state is needed
}


// writes t, theta1, theta2, E, relative and cumulative energy error, once every dt_out
void writeEnergy(double t, double* Y, ofstream& out, int& k, double& err_cum) {

    const double dt_out=0.1; // output interval
    const double eps=1e-12;  // tolerance on the comparison between times

    if(t+eps>=k*dt_out){

        double err = E_err(Y);
        err_cum+=err;
        out << setprecision(6)<<setw(15) << t << setw(15) << Y[0] << setw(15) << Y[2]
            << setw(15) << E(Y) << setw(15) << err << setw(15) << err_cum << endl;
        k++;}
}


void integrate(double* Y, int neq, double tmax, int N, StepFunc stepFunc,
               Callback callback, ofstream& out) {

    double dt = tmax / N;
    double t = 0.0;
    double err_cum=0.0;
    int k=0; // output counter, reset at every integration

    for (int i = 0; i < N; i++) {
        callback(t,Y,out,k,err_cum);   // the callback keeps the integrator independent of what is written
        stepFunc(t, Y, dYdt, dt, neq);
        t += dt;
    }
}

void integrate_linearized(double* Y, int neq, double tmax, int N, StepFunc stepFunc,
                          Callback callback, ofstream& out) {
    double dt = tmax / N;
    double t = 0.0;
    double err_cum = 0.0;
    int k = 0;
    for (int i = 0; i < N; i++) {
        callback(t, Y, out, k, err_cum);
        stepFunc(t, Y, dYdt_linearized, dt, neq);
        t += dt;
    }
}


// evaluates quantity_function on an Ngrid x Ngrid grid of initial conditions
// theta1, theta2 in [0, pi), starting at rest. Output: first row = theta2 values,
// then one row per theta1 (first column = theta1).
void compute_grid(const string& filename, double (*quantity_function)(double*, int), int Ngrid) {

    vector<double> values(Ngrid*Ngrid);

    // the grid points are independent: parallelised with OpenMP if available
    #pragma omp parallel for schedule(dynamic)
    for(int idx=0; idx<Ngrid*Ngrid; idx++){
        int i = idx / Ngrid, j = idx % Ngrid;
        double Y[4];
        init_conditions(Y, M_PI*(static_cast<double>(i)/Ngrid), 0, M_PI*(static_cast<double>(j)/Ngrid), 0);
        values[idx] = quantity_function(Y,4);
    }

    ofstream file(filename);
    file << Ngrid;
    for (int j = 0; j < Ngrid; j++)
        file << setw(15) << M_PI*(static_cast<double>(j)/Ngrid);
    file << endl;

    for(int i=0; i<Ngrid; i++){
        file << setw(15) << M_PI*(static_cast<double>(i)/Ngrid);
        for(int j=0; j<Ngrid; j++)
            file << setw(15) << values[i*Ngrid+j];
        file << endl;
    }
    file.close();
}
