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




void Unif_3D ( double posizione[], double delta, Random *rnd)
{
	double x = posizione [0];
	double y = posizione [1];
	double z = posizione [2];


	posizione [0] = rnd->Rannyu(x-delta , x+delta);

	posizione [1] = rnd->Rannyu(y-delta , y+delta);

	posizione [2] = rnd->Rannyu(z-delta , z+delta);




}
void Gauss_3D ( double posizione[], double delta, Random *rnd)
{



	posizione [0] = rnd->Gauss(posizione[0], delta);

	posizione [1] = rnd->Gauss(posizione[1], delta);

	posizione [2] = rnd->Gauss(posizione[2], delta);




}
int Metropolis (double xn [], FunzioneBase *f, double delta, Random *rnd, bool l) // if i=0 use Unif_3d else Gauss_3d
{
	double k;
	double xp[3];
	xp[0]=xn[0];
	xp[1]=xn[1];
	xp[2]=xn[2];
	if (l==0)
	{
		Unif_3D(xp,delta, rnd);

	}
	else
	{
		Gauss_3D ( xp, delta, rnd );
	}
	k = f->Eval(xp)/f->Eval(xn);
	double alpha = min (1. , k );
	//cout<<xp[0]<<"  "<<xp[1]<<"  "<<xp[2]<<"  "<<xn[0]<<"  "<<f->Eval(xp)<<"   "<<f->Eval(xn)<<"   "<<alpha<<"  "<<k<<endl;

	double r = rnd->Rannyu( 0,1);
	if (r < alpha)
	{

		xn[0]=xp[0];
		xn[1]=xp[1];
		xn[2]=xp[2];

	 return 1;
	}
	else
		return 0;
}

void Simulazione (int M, int N, double x[], double ave[], double av2[], double sum_prog[], double su2_prog[], double err_prog[] )
{


	for (int i=0; i<N; i++)
	{
		ave [i] = 0. ;
		av2 [i] = 0. ;
		sum_prog [i] = 0. ;
		su2_prog [i] = 0. ;
		err_prog [i] = 0. ;


	}

int l=M/N;



for ( int i=0 ; i<N; i++)
	{
 		double sum = 0;
 		int k=0;

		for ( int j=0; j<l; j++)
		{
			k = j+i*l;
			sum = sum + x[k];

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

double devst (double av[], double av2[], int i)
{

	double diff;
	diff = av2[i] - (av[i] * av[i]) ;
	if (i==0)
	        return 0;
	else
		return pow( diff / i , 0.5); // Statistical uncertainty

}
