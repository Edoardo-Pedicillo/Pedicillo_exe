#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.h"


using namespace std;


int main () {

	Random *rnd = new Random ();
	inizRandom ( rnd );

	FunzioneBase *f = new coseno (  );

	int M = 100000; //NUmero di lanci
	int N = 100; //Numero di blocchi
	int l = M/N; //Numero di lancio per blocco


	double r[M]; //contiene i numeri generati casualmente
	double x [N];
	double ave [N];
	double av2 [N];
	double sum_prog [N];
	double su2_prog [N];
	double err_prog [N];
	cout<<"ES.2.1"<<endl<<"Valori settati:"<<endl;
	cout<<"numero lanci: "<<M<<endl<<"Numero blocchi: "<<N<<endl;
	cout<<"Calcolo integrale con distribuzione uniforme "<<endl;

// Integrazione con generatore distribuzione uniforme

	for (int i=0; i<M; i++)
	{
		r[i] = rnd->Rannyu( );

	}
	double appo;



	Simulazione( M, N, f, r, x, ave, av2, sum_prog, su2_prog, err_prog  );

	ofstream fout ("Integrale.dat");

	for (int i=0; i<N; i++)
	{


		fout << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;



	}




	fout.close();
	cout<<"Calcolo integrale con importance sampling "<<endl;
	//Integro usando la distribuzione  -2*(x-1)

	double xmin = 0;
	double xmax = 1;
	double pmax = 1.5;

	FunzioneBase *g = new funzione2( ); // funzione che devo integrare
	ofstream fout3 ("punti.dat");
	for (int i=0; i<M; i++)
	{
		r[i] =rnd->Rannyu();

	}



	Simulazione( M, N, g, r, x, ave, av2, sum_prog, su2_prog, err_prog  );

	ofstream fout2 ("Integrale2.dat");

	for (int i=0; i<N; i++)
	{


		fout2<< x[i]*l <<" "<<sum_prog[i]<<" "<<err_prog[i]<< endl;



	}



	fout2.close();
	rnd->SaveSeed();

	return 0;
}
// SEED: 1642 794 2463 3585
