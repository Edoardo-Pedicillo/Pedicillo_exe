#define _USE_MATH_DEFINES //per usare pi greco

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "classi.h"

using namespace std;

void inizRandom ( Random *); // setta automaticamente i parametri del generatore

double devst (double av[], double av2[], int i);

void Simulazione (int M, int N, FunzioneBase *f, double r[], double x[], double ave[], double av2[], double sum_prog [], double su2_prog[], double err_prog[] ); // del main1.cpp

void Simulazione (int M, int N, FunzioneBase *f, Random *r, double x[], double ave[], double av2[], double sum_prog [], double su2_prog[], double err_prog[] ); // chiede un random r invece di un vett
