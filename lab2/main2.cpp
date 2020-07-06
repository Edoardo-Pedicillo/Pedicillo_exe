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

	FunzioneBase *iden = new identita ();

	int M = 100000;//n cammini
	int N = 100;// N passi
	int l = M/N;
	double **matr = new double* [N]; // i numeri di riga sono i passi e le colonne sono le simulazioni
	double x [N];
	double ave [N];
	double av2 [N];
	double sum_prog [N];
	double su2_prog [N];
	double err_prog [N];

	double posizione [3] = {0,0,0};

	//inizializzazione di matr
	for (int i=0; i<N; i++)
		matr[i]= new double [M];

	cout<<"Esercizio 2.2: "<<endl<<"Parametri simulazione:"<<endl;
	cout<<"numero di simulazioni: "<<M<<endl<<"Numero di passi: "<<N<<endl;

	// caso discreto

	for (int i=0; i<M; i++)
	{
		//posizone iniziale a 0
		posizione[0]=0;
		posizione[1]=0;
		posizione[2]=0;


		for (int j=0; j<N; j++)
		{
			double cas = rnd->Rannyu(0.,3.);

			if (cas > 0 and cas < 1 )
			{
				if (cas > 0.5 )
					posizione [0] ++;

				else
				 	posizione [0] --;
			}
			if (cas > 1 and cas < 2 )
			{
				if (cas > 1.5 )
					posizione [1] ++;

				else
					posizione [1] --;
			}

			if (cas > 2 and cas < 3 )
			{
				if (cas > 2.5 )
					posizione [2] ++;

				else
					posizione [2] --;
			}

			matr[j][i] = posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2] ;



		}
	}





	ofstream fout ("RWDiscreto.dat");

	for (int i=0; i<N; i++)
	{
		//calcolo media
		double sum = 0;

		for ( int j=0; j<M; j++)
		{
				sum += matr[i][j];


		}


		double mean = sum / M;



		Simulazione( M, N, iden, matr[i], x, ave, av2, sum_prog, su2_prog, err_prog  );
		fout << i <<" "<<sqrt(mean)<<" "<<err_prog[N-1]/(2*sqrt(mean))<< endl;



	}


	fout.close();



// Caso continuo

	double teta;
	double phi;


	ofstream fout2 ("RWContinuo.dat");
	for (int i=0; i<M; i++)
	{
		//posizone iniziale a 0
		posizione[0]=0;
		posizione[1]=0;
		posizione[2]=0;


		for (int j=0; j<N; j++)
		{
			teta = acos(1-2*rnd->Rannyu( ));

			phi =  2*M_PI*(rnd->Rannyu());

			posizione [0] += sin(teta) * cos(phi);

			posizione [1] += sin(teta) * sin (phi);

			posizione [2] += cos (teta);



			matr[j][i] = posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2] ;


		}

	}



	for (int i=0; i<N; i++)
	{
		//calcolo media
		double sum = 0;

		for ( int j=0; j<M; j++)
		{
				sum += matr[i][j];

			
		}


		double mean = sum / M;



		Simulazione( M, N, iden, matr[i], x, ave, av2, sum_prog, su2_prog, err_prog  );
		fout2 << i <<" "<<sqrt(mean)<<" "<<err_prog[N-1]/(2*sqrt(mean))<< endl;



	}




	return 0;

}
