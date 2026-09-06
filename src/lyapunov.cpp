#include "simulate.h"
#include "integrators.h"
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

double lyapunov_max_benettin(double* Y, int neq)
{
    const double eps = 1e-8;   // perturbazione iniziale
    const double dt = 0.005;      // passo temporale su cui vado a calcolarmi i singoli esponenti da mediare dopo
    const int Nstep = 50;        // sotto-step RK4 --> h = dt/Nstep = 1e-4
    const int transient = 200;   // iterazioni scartate
    const int Niter = 2000;       //  --> t_medio = 2000*0.005 = 10

    double Y_p[neq]; //vettore perturbato

    for(int i=0;i<neq;i++) Y_p[i] = Y[i];

    // perturbazione iniziale --> mi basta perturbare in una direzione tanto poi si allinea tutto nella direzione di massimo esponente
    Y_p[0] += eps;

    double d0 = eps; //distanza a t=0 delle due soluzioni
    double d, sum = 0.0;

    ofstream nullfile;
    int k = 0;

    //ripeto l esperimento Niter volte, le condizioni iniziali varieranno leggermente ---> faccio poi una media
    for(int n=0; n<Niter; n++){

        // integra entrambe per dt
        integrate(Y , neq, dt, Nstep, RK4Step, dummyCallback, nullfile);
        integrate(Y_p, neq, dt, Nstep, RK4Step, dummyCallback, nullfile);

        // distanza nello spazio delle fasi
        d = 0.0; //inizializzo a 0 tanto conto anche step 0
        for(int i=0;i<neq;i++){
            d += (Y_p[i]-Y[i])*(Y_p[i]-Y[i]);}
        d = sqrt(d);

        // evita crash numerici
        if(d < 1e-14) {
            d = 1e-14;}

        // accumula solo dopo il transiente
        if(n>=transient){
            sum+=log(d/d0);}

        // rinormalizzazione (cuore di Benettin)---> modifico la c.i. successiva in base al passo precedente
        for(int i=0;i<neq;i++){
            Y_p[i] = Y[i] + (Y_p[i]-Y[i])*(d0/d);}

    }

    return sum / ((Niter-transient)*dt);

}









double ftle(double* Y, int neq){
    // esponente di lyapunov a tempo preciso: non usato perche poco preciso

    const double eps = 1e-8;
    const int N = 10000;   // passi totali
    double T=10;
    ofstream nullfile;

    double Y_p[4];

    for(int i=0;i<neq;i++){
        Y_p[i] = Y[i];
    }

    // perturbazione iniziale
    Y_p[0] += eps;

    integrate(Y , neq, T, N, RK4Step, dummyCallback, nullfile);
    integrate(Y_p, neq, T, N, RK4Step, dummyCallback, nullfile);

    double d = 0.0;
    for(int i=0;i<neq;i++)
        d += (Y_p[i]-Y[i])*(Y_p[i]-Y[i]);

    d = sqrt(d);

    return log(d/eps) / T;
}
