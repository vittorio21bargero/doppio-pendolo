#include "integrators.h"
#include "system.h"  // serve per usare RHSFunc
#include <cmath>
#include <string>
#include <vector>

using namespace std;

// definisco i metodi di integrazione
vector<Method> methods = { // gli scruct sono dei tipi quindi al posto di int ... qui avro method
    {"Euler", EulerStep},
    {"RK2", RK2Step},
    {"RK4", RK4Step},
    {"Verlet", position_verlet}
    };


void EulerStep(double t, double *Y, void(*RHSFunc)( double*, double*), double dt, int neq) {
    double rhs[4];
    RHSFunc( Y, rhs);

    for (int k = 0; k < neq; k++) {
        Y[k] += dt * rhs[k]; //per ogni componente ad ogni passo temporale la y varia come il rhs (metodo al primo ordine)
    }
}

//RUNGE-KUTTA 2

void RK2Step (double t, double *Y, void(*RHSFunc)( double*, double*), double h, int neq){

    double Y1[neq], k1[neq], k2[neq];

    RHSFunc(Y,k1); //  k1=R quindi posso direttamente prenderlo dalla funzione

    int i;

    for (i=0;i<neq;i++){
        Y1[i]= Y[i]+0.5*h*k1[i];}

    RHSFunc(Y1,k2); // da capire un po meglio

    for (i=0;i<neq;i++){
        Y[i]+=h*k2[i];}



}



//RUNGE-KUTTA 4

void RK4Step (double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq){

    double Y1[neq], Y2[neq], Y3[neq], k1[neq], k2[neq], k3[neq], k4[neq];
    int i;

    RHSFunc(Y,k1);
    for (i=0;i<neq;i++){
        Y1[i]= Y[i]+0.5*h*k1[i];}

    RHSFunc(Y1,k2);
    for (i=0;i<neq;i++){
        Y2[i]= Y[i]+0.5*h*k2[i];}

    RHSFunc(Y2,k3);
    for (i=0;i<neq;i++){
        Y3[i]= Y[i]+h*k3[i];}

    RHSFunc(Y3,k4);
    for (i=0;i<neq;i++){
        Y[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);}

}

//POSITION VERLET

void position_verlet(double t, double *Y, void(*RHSFunc)(double*, double*), double h, int neq){

    int i;
    double R[neq];

    // calcolo le posizioni a t=h/2
    for (i=0; i<neq; i+=2) //i+=2 per skippare le velocita
        Y[i] += 0.5*h*Y[i+1];

    // calcolo le accelerazioni R(x=x(t=h/2)) di cui poi mi servira solo R[1]
    RHSFunc(Y,R);

    // calcolo velocita e poi posizioni
    for (i=0; i<neq; i+=2) {
        Y[i+1] += h*R[i+1];
        Y[i]   += 0.5*h*Y[i+1];}

    }
