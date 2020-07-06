#include "funzioni.h"

using namespace std;


void inizRandom ( Random *rnd )
{

	int seed[4];
	int p1, p2;

 ifstream Primes("Primes");

   if (Primes.is_open()){

      Primes >> p1 >> p2 ;

	} else cerr << "PROBLEM: Unable to open Primes" << endl;

	Primes.close();

	ifstream input("seed.in");
 	string property;

	if (input.is_open()){

      	while ( !input.eof() ){

        input >> property;

        if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd->SetRandom(seed,p1,p2);
         }
      }

      input.close();

   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


}


double devst (double av[], double av2[], int i)
{

	double diff;
	diff = av2[i] - (av[i] * av[i]) ;
	if (i==0)
	        return 0;
	else
		return pow( diff / i , 0.5); // Statistical uncertainty




}


void Simulazione (int M, int N, FunzioneBase *f, double r[], double x[], double ave[], double av2[], double sum_prog[], double su2_prog[], double err_prog[] )
{

	//settaggio a zero dei vettori

	for (int i=0; i<N; i++)
	{
		ave [i] = 0. ;
		av2 [i] = 0. ;
		sum_prog [i] = 0. ;
		su2_prog [i] = 0. ;
		err_prog [i] = 0. ;
	 	x[i] = i;


	}


	int l = (int) M/N;






	for ( int i=0 ; i<N; i++)
	{
 		double sum = 0;

 		int k=0;



		for ( int j=0; j<l+1; j++)
		{
			k = j+i*l;

			sum = sum + f->Eval(r[k]);


		}

			ave[i] = sum / l;     // r_i
			av2[i] = ave[i] * ave[i];   // r_i^2

		
	}



	for ( int i=0; i<N; i++ )
	{
		for ( int j =0; j<i+1; j++)
		{
			sum_prog[i] += ave[j];  //SUM_{j=0,i} r_j
			su2_prog[i] += av2[j]; // SUM_{j=0,i} (r_j)^2

		}
		sum_prog[i] /= (i+1);  // Cumulative average
		su2_prog[i]/=(i+1); // Cumulative square average


		err_prog[i] = devst( sum_prog , su2_prog, i); // Statistical uncertainty

	}





}



void Simulazione (int M, int N, FunzioneBase *f, Random *r, double x[], double ave[], double av2[], double sum_prog[], double su2_prog[], double err_prog[] )
{

	//settaggio a zero dei vettori

	for (int i=0; i<N; i++)
	{
		ave [i] = 0. ;
		av2 [i] = 0. ;
		sum_prog [i] = 0. ;
		su2_prog [i] = 0. ;
		err_prog [i] = 0. ;
	 	x[i] = i;


	}


	int l=M/N;






	for ( int i=0 ; i<N; i++)
	{
 		double sum = 0;

 		double y;

		for ( int j=0; j<l+1; j++)
		{

			y = f->Eval(r);
			sum = sum + y ;



		}
			ave[i] = sum/l;     // r_i
			av2[i] = ave[i] * ave[i];   // r_i^2


	}



	for ( int i=0; i<N; i++ )
	{
		for ( int j =0; j<i+1; j++)
		{
			sum_prog[i] += ave[j];  //SUM_{j=0,i} r_j
			su2_prog[i] += av2[j]; // SUM_{j=0,i} (r_j)^2

		}
		sum_prog[i] /= (i+1);  // Cumulative average
		su2_prog[i]/=(i+1); // Cumulative square average


		err_prog[i] = devst( sum_prog , su2_prog, i); // Statistical uncertainty

	}





}
