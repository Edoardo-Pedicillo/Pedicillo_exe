#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.h"

//Includendo "classi.h" mi da codici di errore

using namespace std;



int main () {

	Random *rnd = new Random ();
	inizRandom ( rnd );



	int M = 10000;
	int N [4] = {1,2,10,100};


	// Dado normale

	ofstream fout1 ("dadoUnif.dat");

	FunzioneBase *f = new unifDistribution(  );

	for (int i=0 ; i<4; i++)
		{
			Dadi ( fout1, f, rnd,  M, N[i] ); // simula i lanci
		}

fout1.close();
	//Dado esponenziale

	ofstream fout2 ("dadoExp.dat");

	double lambda=1;

	FunzioneBase *esp = new expDistribution ( lambda );

	for (int i=0 ; i<4; i++)
		{
			Dadi ( fout2, esp, rnd,  M, N[i] ); // simula i lanci
		}



	fout2.close ( );

	//Dado Lorentziano

	ofstream fout3("dadoLor.dat");

	double gamma=1;

	FunzioneBase *lor = new LorDistribution ( gamma );

	for (int i=0 ; i<4; i++)
		{
			Dadi ( fout3 , lor , rnd,  M, N[i]); // "Dadi" calcola Sn e ne restituisce il valore ad un file 
		}



	fout3.close ( );


	return 0;
}
