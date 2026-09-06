#include "simulate.h"
#include "energy.h"
#include "system.h"
#include "integrators.h"
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;



void init_conditions(double* Y, double th1, double w1, double th2, double w2) {
    Y[0] = th1;
    Y[1] = w1;
    Y[2] = th2;
    Y[3] = w2;
}


void dummyCallback(double t, double* Y, ofstream& out, int& k,double& err_cum) {
    // non fa nulla
}




// funzione callback che scrive dati su file
void writeData(double t, double* Y,ofstream& out) {
    out << setprecision(4)<<setw(10) << t << setw(10) << Y[0] << setw(10) << Y[2]
          << setw(10) << E(Y) <<  setw(13) << E_err(Y) << endl;
}


void writeData1(double t, double* Y,ofstream& out,int& k,double& err_cum) {

    const double dt_out=0.1; //voglio stampare solo ogni dt=0.1
    const double eps=1e-12;  //tolleranza

    if(t+eps>=k*dt_out){

        double err = E_err(Y);
        err_cum+=E_err(Y);
        out << setprecision(4)<<setw(15) << t << setw(15) << Y[0] << setw(15) << Y[2]
          << setw(15) << E(Y) <<  setw(15) << err << setw(15) << err_cum << endl;
          k++;}
}



void integrate(double* Y, int neq, double tmax, int N,
               void(*stepFunc)(double, double*, void(*)(double*, double*),double,int),
               void(*callback)(double, double*,ofstream&,int& k,double& err_cum),
               ofstream& out) {

    double dt = tmax / N;
    double t = 0.0;
    double err_cum=0.0;


    int k=0; //ad ogni ciclo di integrazione voglio che si annulli k

    for (int i = 0; i < N; i++) {
        callback(t,Y,out,k,err_cum);          // ad ogni step scrivo i dati con una funzione di callback: molto modulare
        stepFunc(t, Y, dYdt, dt, neq);
        t += dt;
    }
}

void integrate_linearized(double* Y, int neq, double tmax, int N,
               void(*stepFunc)(double, double*, void(*)(double*, double*), double, int),
               void(*callback)(double, double*, ofstream&, int&, double&),
               ofstream& out) {
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


//string  mi permette di creare direttamente il file all'interno della funzione
void compute_grid(const string& filename,
                   double (*quantity_function)(double*,int))
                    {
    ofstream file(filename);
    double Y[4];
    int Ngrid = 100;

    // riga di header: valori di theta2 (colonne)
    file << Ngrid;
    for (int j = 0; j < Ngrid; j++)
        file << setw(15) << M_PI*(static_cast<double>(j)/Ngrid);
    file << endl;

    for(int i=0; i<Ngrid; i++){
        double th1 = M_PI*(static_cast<double>(i)/Ngrid);
        file << setw(15) << th1;   // prima colonna: valore di theta1 (riga)
        for(int j=0; j<Ngrid; j++){
            double th2 = M_PI*(static_cast<double>(j)/Ngrid);
            init_conditions(Y, th1, 0, th2, 0);
            double value = quantity_function(Y,4);
            file << setw(15) << value;
        }
        file << endl;
    }
    file.close();
}
