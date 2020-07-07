#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

class individuo {

public:

individuo(int N_cities_m, Random mrnd, double **posizioni_m);
bool check( ); // 0: non va bene , 1: ok
void mutation( double t);
double Pbc(double r);
double L2 ( );
double Get (int n) {return ind[n];};
void Set (int n, int city) {ind [n]=city;}
void Swap (int a, int b) {std::swap(ind[a],ind[b]);};
void copy(individuo i2);
void crossover( individuo i2);
bool operator< ( individuo individuo1){return L2() < individuo1.L2();}
private:
  vector <int> ind;
  int N_cities;
  Random rnd;
  double **posizioni;
  double prob_mutation=0.01;
  double prob_crossover=0.7;
  int p=50;
  int n=5; // mutation (2) : di quanto traslare
  int m2=2; // quante città traslare (m<N_cities-1)
  int m3 = 10; // mutation[3] m<N/2


};


class population {
public:
  population ( int Npop_m , int  N_cities_m , Random mrnd, double **posizioni_m);
  void start();
  void totalcheck();
  void order ();
  bool converge(); // 0:si e 1:no
  int selection();
  double L2 (int i) {return pop[i].L2();};
  double Pbc(double r);
  int Get(int i,int j) {return pop[j].Get(i);};
  void crossover(int genitore1, int genitore2);
  void genetic (double t);
  double average();



private:
  int Npop, N_cities;
  double **posizioni;
  double prob_mutation=0.10;
  double prob_crossover=0.70;
  Random rnd;
  vector <individuo> pop  ;
  vector <individuo> new_pop;
  int p=50;
  int n=5; // mutation (2) : di quanto traslare
  int m2=6; // quante città traslare (m<N_cities-1)
  int m3 = 10; // mutation[3] m<N/2
};
