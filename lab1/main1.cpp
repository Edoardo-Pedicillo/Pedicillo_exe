#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.h"


using namespace std;


int main () {


	int M = 10000; //NUmero di lanci
	int N = 100; //Numero di blocchi
	int l = M/N; //Numero di lancio per blocco


	double r[M]; //contiene i numeri generati casualmente

	Random *rnd = new Random ();

	//Settaggio del generatore Random

	inizRandom( rnd );


	double x [N];
	double ave [N];
	double av2 [N];
	double sum_prog [N];
	double su2_prog [N];
	double err_prog [N];




	//Inizio simulazione es 1.1.1

	for (int i=0; i<M; i++)
	{
		r[i] = rnd->Rannyu( );

	}

	FunzioneBase *f = new unifDistribution( );

	Simulazione( M, N, f, r, x, ave, av2, sum_prog, su2_prog, err_prog  );



	// Passo i dati al file "data1.dat"

	ofstream fout ("data1.dat");

	for (int i=0; i<N; i++)
	{


		fout << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;



	}




	fout.close();


	// Inizio simulazione es 1.1.2

	FunzioneBase* deviazione = new dev();


	Simulazione ( M, N, deviazione, r , x , ave , av2 , sum_prog ,su2_prog , err_prog  );




	ofstream fout1 ("data2.dat");

	for (int i=0; i<N; i++)
	{


		fout1 << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;



	}




	fout1.close();



	//Esercizio 1.1.3

	M = 100 ; //sottointervalli
	int nmis=10000; //lanci ad ogni ciclo
	double E ; //valore atteso
	double subint[M]; //Tiene conto dei numeri generati che cadono in ciascuno degli M sottointervalli
	double ciclo [M];
	double csi [M];
	double n;

	// setto gli array

	for (int i=0; i<M+1; i++)
	{

		subint[i] = 0;
		ciclo[i] = i;
		csi[i] = 0;
	}



	for (int time=1; time < M+1; time++)
	{
		for (int i=0; i<M+1; i++)
		{

			subint[i] = 0;

		}

		for (int i=0; i< nmis+1 ; i++) // fa partire il generatore 10^4 volte
		{


			n = rnd->Rannyu( );


			for (int j=1; j<M+1; j++) // Colloca ciascuno dei numeri generati nell'intervallo a cui appartengono
			{
				if (n < (double) j/M and n > (double)( j - 1 )/ M)
				{
					subint[j-1] = subint[j-1]+1 ;
				}


			}




		}




		E = (double) nmis/M;



		for (int j=0; j+1<M; j++)
		{
			csi[time] += chi2 (subint , E, j); //la funzione chi2 calcola solo l'addendo della formula del chi-quadro

		}


	}
	ofstream fout2 ("data3.dat");

	for (int j=1; j<M+1; j++)
		{
			fout2<< ciclo[j] << " " << csi[j] <<endl;

		}
	fout2.close( );

	return 0;
}
