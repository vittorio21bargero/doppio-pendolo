// Double pendulum: comparison of integrators and characterisation of chaos.
//
// Usage:  ./double_pendulum [experiment]
//   energy       relative energy error vs time, for each method and time step
//   convergence  global error vs step size on the linearized system (exact solution known)
//   lyapunov     largest Lyapunov exponent, FTLE and flip time vs the initial angle
//   grids        flip number and Lyapunov exponent on a 100x100 grid of initial conditions
//   all          everything above (default)
// Results are written as text files in data/; scripts/plot_results.py turns them into figures.

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <filesystem>
#include "integrators.h"
#include "system.h"
#include "energy.h"
#include "simulate.h"
#include "lyapunov.h"
#include "flip_time.h"

using namespace std;

const int neq = 4;


// ENERGY ERROR vs t, for the 4 methods and dt = 1e-2 ... 1e-6
// initial condition: both arms horizontal, at rest
void energy_experiment() {

    filesystem::create_directories("data/energy");
    double Y[neq];
    double tmax = 10.0;

    for (int i = 2; i <= 6; i++) {
        double dt = pow(10.0, -i);
        int N = static_cast<int>(round(tmax / dt));

        for (const Method& m : methods) {
            ofstream out("data/energy/" + m.name + "_dt_1e-" + to_string(i) + ".txt");
            out << "# method: " << m.name << "   dt = " << dt << endl;
            out << "#      t         theta1         theta2              E             dE         dE_cum" << endl;

            init_conditions(Y, 0.5*M_PI, 0.0, 0.5*M_PI, 0.0);
            E_0 = E(Y);
            integrate(Y, neq, tmax, N, m.func, writeEnergy, out);
        }
        cout << "energy: dt = " << dt << " done" << endl;
    }

    // ENERGY ERROR vs t for increasing amplitude theta1 = theta2 = pi*i/M, RK4, dt = 1e-5
    int M = 6;
    for (int i = 1; i < M; i++) {
        ofstream out("data/energy/RK4_amplitude_" + to_string(i) + "of" + to_string(M) + ".txt");
        double theta = M_PI*(static_cast<double>(i)/M);
        out << "# RK4, dt = 1e-5, theta1 = theta2 = " << theta << endl;
        init_conditions(Y, theta, 0, theta, 0);
        E_0 = E(Y);
        integrate(Y, neq, tmax, 1000000, RK4Step, writeEnergy, out);
    }
    cout << "energy: amplitude scan done" << endl;
}


// CONVERGENCE STUDY: on the linearized system theta1(t) is known exactly,
// so the global error at t = 10 can be measured as a function of h = 2^-n
void convergence_experiment() {

    filesystem::create_directories("data/convergence");

    double tmax = 10;
    double g = 9.81;
    double L = 1;
    double theta0 = 1e-1;
    double theta1_exact = theta1_theoretical(tmax, theta0, g, L);

    ofstream nullfile; // unused stream for the empty callback

    for (const Method& m : methods) {
        ofstream out("data/convergence/" + m.name + ".txt");
        out << "# method: " << m.name << endl;
        out << "#           h      theta1_num    theta1_exact       abs_error" << endl;

        for (int n = 2; n <= 20; ++n) {
            double h = pow(2.0, -n);
            int N = static_cast<int>(round(tmax / h));

            double Y[neq];
            init_conditions(Y, theta0, 0.0, theta0, 0.0);
            integrate_linearized(Y, neq, tmax, N, m.func, dummyCallback, nullfile);

            double err = fabs(Y[0] - theta1_exact);
            out << setprecision(10) << setw(16) << h << setw(16) << Y[0]
                << setw(16) << theta1_exact << setw(16) << err << endl;
        }
        cout << "convergence: " << m.name << " done" << endl;
    }
}


// LYAPUNOV EXPONENT, FTLE and FLIP TIME along the diagonal theta1 = theta2 = theta in [0, pi]
void lyapunov_experiment() {

    filesystem::create_directories("data/lyapunov");

    const int Nstep = 100;
    vector<double> lam_ftle(Nstep+1), lam_max(Nstep+1), t_flip(Nstep+1);

    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j <= Nstep; j++) {
        double theta = M_PI * (static_cast<double>(j)/Nstep);
        double Y[neq];

        init_conditions(Y, theta, 0, theta, 0);
        lam_ftle[j] = ftle(Y, neq);

        init_conditions(Y, theta, 0, theta, 0);
        lam_max[j] = lyapunov_max_benettin(Y, neq);

        init_conditions(Y, theta, 0, theta, 0);
        t_flip[j] = flip_time(Y, neq);
    }

    ofstream out("data/lyapunov/lyapunov_vs_theta.txt");
    out << "#   theta          FTLE    lambda_max     flip_time" << endl;
    for (int j = 0; j <= Nstep; j++)
        out << setw(10) << M_PI*(static_cast<double>(j)/Nstep) << setw(14) << lam_ftle[j]
            << setw(14) << lam_max[j] << setw(14) << t_flip[j] << endl;
    cout << "lyapunov: done" << endl;
}


double flip_number_d(double* Y, int neq) { return flip_number(Y, neq); } // adapter for compute_grid

// MAPS over the plane of initial angles (theta1, theta2), both arms at rest
void grids_experiment() {

    filesystem::create_directories("data/grids");

    compute_grid("data/grids/flip_number.txt", flip_number_d);         // T = 50
    cout << "grids: flip number done" << endl;

    compute_grid("data/grids/lyapunov.txt", lyapunov_max_benettin);   // T = 10, Benettin
    cout << "grids: lyapunov done" << endl;
}


int main(int argc, char** argv) {

    string which = (argc > 1) ? argv[1] : "all";

    if (which == "energy"      || which == "all") energy_experiment();
    if (which == "convergence" || which == "all") convergence_experiment();
    if (which == "lyapunov"    || which == "all") lyapunov_experiment();
    if (which == "grids"       || which == "all") grids_experiment();

    if (which != "energy" && which != "convergence" && which != "lyapunov" &&
        which != "grids" && which != "all") {
        cerr << "unknown experiment '" << which
             << "'. Use: energy | convergence | lyapunov | grids | all" << endl;
        return 1;
    }
    return 0;
}
