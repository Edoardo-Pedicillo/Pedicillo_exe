#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.h"

//Includendo "classi.h" mi da codici di errore

using namespace std;


int main ()
{
	
	int nRighe = 100000;
	double L = 0.5; // lunghezza ago
	double d = 1.0;	//distanza tra righe

	double xmax = (( double ) nRighe - 1.0) * d ; // Supponiamo che "il pavimento " è quadrato quindi ymax = xmax
	int M = 50000;
	int N = 100;
	int l = M/N;
	Random *rnd = new Random ();

	inizRandom ( rnd );

	double x [N];
	double ave [N];
	double av2 [N];
	double sum_prog [N];
	double su2_prog [N];
	double err_prog [N];


	FunzioneBase *f = new Buffon ( nRighe, L, d, xmax, M, N );

	Simulazione (M, N, f, rnd,  x, ave, av2, sum_prog , su2_prog, err_prog );



	ofstream fout ("Buffon.dat");
	ofstream fout0 ("Buffon0.dat");

	for (int i=0; i<N; i++)
	{


		fout0 << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;



	}

	for (int i=0; i<N; i++)
	{

		sum_prog[i] = 1 / sum_prog [i]; // In "Simulazione" calcolo 1/pi per poter riciclare il codice già scritto per es 1.1.1

		err_prog [i] = err_prog[i] * (sum_prog[i] * sum_prog[i]); // uso la propagazione degli errori per calcolare l'errore sulla misura di pi sapendo l'errore sulla misura di 1/pi

		fout << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;



	}




	fout.close();
	fout0.close();

	return 0;
}
