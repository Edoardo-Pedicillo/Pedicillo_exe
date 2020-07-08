#include "random.h"

//random
int seed[4];
Random rnd;

double hbar = 1;
double m = 1;
double dt = 0.0005;

int Nthrow = 1000000;
int Nblocks = 100;



void Input();
double psi (double x, double sigma, double mu);
double psi2 (double x, double sigma, double mu);
double pot (double x);
double energy (double x, double sigma, double mu);
double Unif_3D ( double posizione);
int Metropolis (double &x,double delta, double sigma, double mu);
void Simulazione (int M, int N, double x[], double ave[], double av2[], double sum_prog[], double su2_prog[], double err_prog[] );
double devst (double av[], double av2[], int i);
double Minimize (double parametri[] , double delta, double delta1, double varE, double t);
