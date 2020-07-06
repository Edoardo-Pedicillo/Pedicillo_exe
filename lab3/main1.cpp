#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "funzioni.h"


using namespace std;


int main () {

// Parametri prezzi
	double S_0 = 	100;
	double T = 1;
	double K = 100;
	double r = 0.1;
	double sigma = 0.25;

// Parametri esperimento

	int M = 10000; //NUmero di lanci
	int N = 100; //Numero di blocchi
	int l = M/N; //Numero di lancio per blocco

	Random *rnd = new Random ();
	inizRandom( rnd );

double x [N];
double ave [N];
double av2 [N];
double sum_prog [N];
double su2_prog [N];
double err_prog [N];

cout<<"_____ES.4_____"<<endl<<endl<<":::::Parametri:::::"<<endl;
cout<<endl<<"S_0: "<<std::setw(10)<<S_0<<endl<<"T: "<<setw(10)<<T<<endl<<"K: "<<setw(10)<<K<<endl<<"r: "<<setw(10)<<r<<endl<<"s: "<<setw(10)<<sigma<<endl;
cout<<"M: "<<setw(10)<<M<<endl<<"N: "<<setw(10)<<N<<endl;
// Part 1

FunzioneBase *f = new PriceCall( S_0 , T , K , r , sigma, 1. );
Simulazione( M, N, f, rnd, x, ave, av2, sum_prog, su2_prog, err_prog  );

ofstream fout ("S1.dat");

for(int i=0; i<N; i++)
{
		fout << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;
}
fout.close();


//Part 2

	int stepT = 100;

	FunzioneBase *g = new PriceCall2( S_0 , T , K , r , sigma, stepT, 1. );
	Simulazione( M, N, g, rnd, x, ave, av2, sum_prog, su2_prog, err_prog  );

	ofstream fout2 ("S2.dat");

	for(int i=0; i<N; i++)
	{
			fout2 << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;
	}
	fout2.close();

	// Part 3

	FunzioneBase *f1 = new PriceCall( S_0 , T , K , r , sigma, -1. );
	Simulazione( M, N, f1 , rnd, x, ave, av2, sum_prog, su2_prog, err_prog  );

	ofstream fout3 ("P1.dat");

	for(int i=0; i<N; i++)
	{
			fout3 << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;
	}
	fout3.close();


	//Part 4



		FunzioneBase *g1 = new PriceCall2( S_0 , T , K , r , sigma, stepT, -1. );
		Simulazione( M, N, g1 , rnd, x, ave, av2, sum_prog, su2_prog, err_prog  );

		ofstream fout4 ("P2.dat");

		for(int i=0; i<N; i++)
		{
				fout4 << x[i]*l <<" "<<sum_prog[i] <<" "<<err_prog[i] << endl;
		}
		fout4.close();
	return 0;
}
