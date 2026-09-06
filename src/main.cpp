#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "integrators.h"
#include "system.h"
#include "energy.h"
#include "simulate.h"
#include "lyapunov.h"
#include "flip_time.h"
#include <filesystem>

using namespace std;





int main() {

    vector<Method> methods = {
    {"Euler",  EulerStep},
    {"RK2",    RK2Step},
    {"RK4",    RK4Step},
    {"Verlet", position_verlet}
};
    const int neq = 4;
    double Y[neq];

    double tmax = 10.0;


    //ENERGY_ERROR(t,dt) per i 4 metodi

    //energia iniziale
    init_conditions(Y, 0.5*M_PI, 0, 0.5*M_PI, 0);
    E_0 = E(Y);


    for (int i = 2; i <= 6; i++) {

        double dt = pow(10.0, -i);
        int N = static_cast<int>(tmax / dt); //esplicito conversione double ---> int

        string filename = "data/error(dt)/dati_dt_" + to_string(i) + ".txt"; // mi serve string type per concatenare il tutto
        ofstream angles_energy_file(filename); //apre direttamente il file senza bisogno di fdata.open()
        angles_energy_file << "#        dt = " << dt << endl<<endl;

        // ciclo sui metodi di integrazione
        for (size_t j = 0; j < methods.size(); ++j) {

            Method& m = methods[j];


            // header utile per gnuplot---> togliere
            angles_energy_file<<endl;
            angles_energy_file << "#        Metodo: " << m.name << endl<<endl;

            angles_energy_file<< "#       t     theta1     theta2       E          dE" << endl;
            angles_energy_file << endl;

            init_conditions(Y, 0.5*M_PI, 0.0, 0.5*M_PI, 0.0);
            integrate(Y, neq, tmax, N, m.func,writeData1, angles_energy_file);
    }

    angles_energy_file.close();
}







   //ENERGY_ERROR(t,condizione iniziale) , dt=10^-4  , metodo = RK4

   int M=6;
   for(int i=1;i<M;i++){

        string filename = "data/ENERGY/ENERGY_init_cond/ENERGY_init_cond_dt_-3_tmax_1000_" + to_string(i) + ".txt"; // mi serve string type per concatenare il tutto
        ofstream angles_energy_file_2(filename); //apre direttamente il file senza bisogno di fdata.open()
        init_conditions(Y,M_PI*(static_cast<double>(i)/M),0,M_PI*(static_cast<double>(i)/M),0);
        E_0=E(Y);
        integrate(Y,neq,tmax,1e6,RK4Step,writeData1,angles_energy_file_2);
    }*/





    //ENERGY_ERROR(t,condizione_iniziale)





    //LYAPUNOV +FLIP TIME
    init_conditions(Y,M_PI/2,0,M_PI/2,0);
    compute_grid("data/lyapunov.txt", ftle);
    init_conditions(Y,M_PI/2,0,M_PI/2,0);
    compute_grid("data/fliptime.txt", flip_time);*/




int Nstep = 100;                 // numero di condizioni iniziali tra 0 e pi



ofstream out("data/LYAPUNOV/lyapunov_theta_tmax_10.txt");


out << "# theta   FTLE   Lyapunov_max\n";

for(int j=0; j<=Nstep; j++){
    double theta = M_PI * (static_cast<double>(j)/Nstep);   // condizione iniziale simmetrica

    // --- FTLE
    init_conditions(Y, theta, 0, theta, 0);
    double lambda_ftle = ftle(Y, neq);

    // --- Lyapunov (Benettin)
    init_conditions(Y, theta, 0, theta, 0);
    double lambda_max = lyapunov_max_benettin(Y, neq);

    // --- flip time
    init_conditions(Y, theta, 0, theta, 0);
    double t_flip = flip_time(Y, neq);

    out << theta << "   " << lambda_ftle << "   "
        << lambda_max << "   " << t_flip << "\n";

    cout << "theta = " << theta
         //<< "  FTLE = " << lambda_ftle
         //<< "  lyap = " << lambda_max << endl;   // per seguire l'avanzamento
}
out.close();*/





// STUDIO DI CONVERGENZA



    ofstream out("data/ERROR_CONVERGENCE/error_convergence.txt");

    double tmax = 10;
    double g  = 9.81;
    double L = 1;
    double theta0 = 1e-1;
    double theta1_exact = theta1_theorical(tmax, theta0, g, L);


    ofstream nullfile; // stream inutilizzato per il callback vuoto



    // ciclo sui metodi di integrazione
    for(size_t j = 0; j < methods.size(); ++j){
        out << "# metodo: " << methods[j].name << "\n";
        out << "# h    theta1_exp    theta1_teo    errore_abs\n";

        // ciclo sugli step di integrazione dt
        for(int n = 2; n <= 20; ++n){
            double h = pow(2.0, -n);
            int N = static_cast<int>(round(tmax / h)); // numero di passi

            // condizioni iniziali: theta1=theta2=theta0, velocita nulle
            double Y[4];

            init_conditions(Y, theta0, 0.0, theta0, 0.0); // adatta all'ordine del TUO stato!

            // integrazione fino a tmax, senza scrivere nulla durante il moto
            integrate_linearized(Y, neq, tmax, N, methods[j].func , dummyCallback, nullfile);

            double theta1_num = Y[0];  // <-- theta1 finale: adatta l'indice se serve
            double err = fabs(theta1_num - theta1_exact);

            out << h << "   " << theta1_num << "   "
                << theta1_exact << "   " << err << "\n";
        }
        out << "\n\n"; // due righe vuote = nuovo blocco index per gnuplot
    }
    out.close();*/



// NUMERO DI FLIP (T = 50)


compute_grid("data/LYAPUNOV_FLIPTIME/fliptime_grid.txt", flip_number);


// LYAPUNOV (T=10,benettin)

compute_grid("data/LYAPUNOV_FLIPTIME/lyapunov_grid.txt", lyapunov_max_benettin);


return 0;
}
