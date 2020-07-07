#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"
#include "classi.h"
using namespace std;

void Init();
int seed[4];
Random rnd;
double pi = 3.1415926535;

int main(){


  int Ncities=32;
  int Npop=100;
  double T = 10;
  double tassoT=0.9999;
  ofstream fout ("positions_c.dat");
  ofstream fout1 ("solution_c.dat");
  ofstream fout2("best_c.dat");
  ofstream fout3("average_c.dat");
  Init();
  double **posizioni= new double *[Ncities];

  //circle
  cout<<"______ES 10_____"<<endl;
  cout<<endl<<"Parte 1"<<endl<<"Parametri:"<<endl;
  cout<<endl<<"Numero città: "<<Ncities<<endl<<"Numero individui: "<<Npop<<endl;
  cout<<"Temperatura: "<<T<<endl<<"Tasso: "<<tassoT<<endl;
  cout<<endl<<"CIRCLE"<<endl;
  double theta=0;
  for (int i=0; i<Ncities; i++)
  {
      posizioni[i]=new double [2];
      theta =rnd.Rannyu(0,2*pi); //rnd.Rannyu(theta,(theta+2*pi)/2);

      posizioni[i][0] = cos(theta);
      posizioni [i][1] = sin(theta);
      fout<<posizioni[i][0] <<" "<<posizioni[i][1] <<endl;
  }
  population p (Npop,Ncities,rnd,posizioni);
  p.start();


  int k=0;
  do
  {
    p.genetic(T);
    fout2<<k<<" "<<p.L2(0)<<endl;
    fout3<<k<<" "<<p.average()<<endl;
    k++;
    T *= tassoT;
  }while(k<5000);/// Faccio girare più volte il programma e abbasso il limite


  for (int i=0; i<Ncities;i++)
  {
    fout1<<posizioni[p.Get(i,0)-1][0]<<" "<<posizioni[p.Get(i,0)-1][1]<<endl;
  }

  fout.close();
  fout1.close();
  fout2.close();
  fout3.close();

  //Square
  ofstream fouts ("positions_s.dat");
  ofstream fout1s ("solution_s.dat");
  ofstream fout2s("best_s.dat");
  ofstream fout3s("average_s.dat");

  cout<<"SQUARE"<<endl;


  for (int i=0; i<Ncities; i++)
  {
      posizioni[i]=new double [2];
      posizioni[i][0] = rnd.Rannyu();
      posizioni [i][1] = rnd.Rannyu();
      fouts<<posizioni[i][0] <<" "<<posizioni[i][1] <<endl;
  }

  population ps (Npop,Ncities,rnd,posizioni);
  ps.start();

  k=0;
  T=150000;
  do
  {

    ps.genetic(T);
    fout2s<<k<<" "<<ps.L2(0)<<endl;
    fout3s<<k<<" "<<ps.average()<<endl;
    k++;
    T *= tassoT;
  }while(k<30000);





  for (int i=0; i<Ncities;i++)
  {

    fout1s<<posizioni[ps.Get(i,0)-1][0]<<" "<<posizioni[ps.Get(i,0)-1][1]<<endl;
  }



  fouts.close();
  fout1s.close();
  fout2s.close();
  fout3s.close();
  rnd.SaveSeed();
  rnd.SaveSeed();
  return 0;
}

void Init()
{

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
            rnd.SetRandom(seed,p1,p2);
         }
      }

      input.close();

   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


}
