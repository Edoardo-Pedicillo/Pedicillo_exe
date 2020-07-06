#define _USE_MATH_DEFINES //per usare pi greco

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "classi.h"

using namespace std;

void inizRandom ( Random *); // setta automaticamente i parametri del generatore

void Unif_3D ( double posizione[], double delta, Random *rnd);

void Gauss_3D ( double posizione[], double delta, Random *rnd);

int Metropolis (double xn [], FunzioneBase *f, double delta, Random *rnd, bool l);

void Simulazione (int M, int N, double x[], double ave[], double av2[], double sum_prog[], double su2_prog[], double err_prog[] );

double devst (double av[], double av2[], int i);
