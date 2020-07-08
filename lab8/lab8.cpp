#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "lab8.h"


using namespace std;

int main () {

  double delta=0.005;
  double parametri [2]= {0.630825,0.806631}; //sigma e mu
  double delta1=0.01;
  double varE = 0.00037;
  double x=0;
  int acc=0;

  double ave [Nblocks];
  double av2 [Nblocks];
  double sum_prog [Nblocks];
  double su2_prog [Nblocks];
  double err_prog [Nblocks];
  double En[Nthrow];

  cout<<"____Monte Carlo Metodo Variazionale_____"<<endl;
  cout<<endl<<"delta1: "<<delta<<endl<<"delta2: "<<delta1<<endl<<"sigma: "<<parametri[0]<<endl<<"mu: "<<parametri[1]<<endl;
  cout<<"dE: "<<varE<<endl<<"dt: "<<dt<<endl;
  ofstream fout ("energy.dat");
  ofstream fout1("psi.dat");
  Input();
  double j= Minimize(parametri, delta, delta1, varE,0.1);

  cout<<"sigma :"<<parametri[0]<<endl<<"mu: "<<parametri[1]<<endl;
  delta=2.5;

  // Istogramma
  for (int i=0; i<Nthrow; i++)
  {

  		acc += Metropolis (x, delta, parametri[0], parametri[1]);
      En[i] = energy(x,parametri[0], parametri[1]);
      fout1<<x<<endl;
  }
  cout<<"acceptance: "<<(double)acc/Nthrow<<endl;

  Simulazione ( Nthrow,Nblocks, En, ave, av2, sum_prog, su2_prog, err_prog );


  // grafico energia
  for (int i=1; i<Nblocks+1; i++)
  {


  	fout << i <<" "<<sum_prog[i-1] <<" "<<err_prog[i-1] << endl;



  }
  fout.close();
  fout1.close();

  return 0;
}


void Input()
{
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();
}

double psi (double x,double sigma,double mu)
{
  return exp(- (x-mu)*(x-mu)/(2*sigma*sigma))+ exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double psi2 (double x,double sigma, double mu)
{
  double y1 = exp (-((mu+x)*(mu+x)/(2*sigma*sigma)))*(mu*mu-sigma*sigma+x*x+2*mu*x)/pow(sigma,4);
  double y2 = exp (-((mu-x)*(mu-x)/(2*sigma*sigma)))*(mu*mu-sigma*sigma+x*x-2*mu*x)/pow(sigma,4);

  return y1+y2;
}

double pot (double x)
{
  return pow(x,4)- 5*pow(x,2)/2;
}

double energy (double x, double sigma, double mu)
{
  return (- (hbar*hbar/(2*m))*psi2(x,sigma,mu)+pot(x)*psi(x,sigma,mu))/psi(x,sigma,mu);
}

double Unif_3D ( double posizione, double d)
{
	return rnd.Rannyu(posizione-d , posizione+d);

}

int Metropolis (double &x,double delta, double sigma, double mu)
{
	double k;


	double x1 = Unif_3D(x,delta);



	k = pow( psi (x1, sigma, mu )/ psi (x,sigma, mu) , 2);

	double alpha = min (1. , k );

	double r = rnd.Rannyu( 0,1);

	if (r < alpha)
	{

		x = x1;

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
		return pow( diff / (i-1 ), 0.5); // Statistical uncertainty

}

double Minimize (double parametri[] , double delta, double delta1 , double varE, double t)
{
  double sigma = parametri[0];
  double mu = parametri[1];
  double sigma1 = Unif_3D(sigma,delta1);
  double mu1 = Unif_3D(mu,delta1);
  double en=0, en1=0;
  double beta = 1./t;
  double p;

  double x=0;
  int acc =0;
  //calcolo energia della funzioe d'onda iniziale
  for (int i=0; i<Nthrow; i++)
  {

  	acc += Metropolis (x, delta, sigma, mu);
    en +=energy(x,sigma,mu);

  }

  en = en/Nthrow;
  //calcolo energia della funzione d'onda proposta
  x=0;
  acc=0;

  for (int i=0; i<Nthrow; i++)
  {

  	acc += Metropolis (x, delta, sigma1, mu1);
    en1 += energy(x,sigma1,mu1);

  }

  en1 = en1/Nthrow;

  // Selezione

  p =  exp(- beta*(en1-en));

  if (en1<en and rnd.Rannyu()<p )
  {



    if(fabs(en-en1)>varE)
    {
      parametri[0]=sigma1;
      parametri[1]=mu1;
      if(t-dt<0)
        cout<<"Errore temperatura < 0"<<endl;
      Minimize(parametri,delta,delta1,varE,t-dt);
    }

    else
      return 0;

  }
  else
  {

    Minimize(parametri,delta,delta1,varE,t-dt);

  }

  return 0;
}
