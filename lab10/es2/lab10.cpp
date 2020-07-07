#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "random.h"
#include "classi.h"
#include "mpi.h"
using namespace std;

void Init(int rank);
int seed[4];
Random rnd;
double pi = 3.1415926535;

int main(int argc, char* argv[])
{
  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat1;
  MPI_Request req1;

  int Ncities=32;
  int Npop=100;

  double T = 50;
  double tassoT=0.9999;

  int send [Ncities];
  int recive [Ncities];
  string s = std::to_string(rank);

  ofstream fout ("positions_c"+s+".dat");
  ofstream fout1 ("solution_c"+s+".dat");
  ofstream fout2("best_c."+s+".dat");
  ofstream fout3("average_c"+s+".dat");
  ifstream pos ("positions_s.dat");

  Init(rank);
  double **posizioni= new double *[Ncities];

  // Acquisizione dati in "positions_s.dat" dove sono contenute le posizioni delle città

  for (int i=0 ; i<Ncities; i++)
  {
    posizioni[i]=new double [2];
    pos>>posizioni[i][0]>>posizioni[i][1];
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



    // Copia della lista delle città della soluzione migliore in send
    for(int i=0;i<Ncities;i++)
    {
      send[i]=p.Get(i,0);
    }

    int destinatario [4] ;


    // viene riempito (dal nodo 0) il vettore "destinatario" di dimensione 4 che contiene i nodi
    // opportunamente rimescolati: questo ci servirà per andare a far scambiare le soluzioni migliori
    if (rank==0)
    {

      for (int i=0; i<size; i++)
      {
        destinatario[i]=i;
      }

      for(int i=0; i<size; i++)
      {
        int x = (int)rnd.Rannyu(0,size);
        std::swap(destinatario[i],destinatario[x]);
      }


    }

    // passo il vettore "destinatario" a tutti gli altri nodi

    MPI_Bcast(&destinatario,4,MPI_INT,0,MPI_COMM_WORLD);


    // ogni 100 passi le prime due coppie e le ultime due in "destinatario" si scambiano la soluzione migliore
    if(k%100==0)
    {
      if(rank==destinatario[0])
      {

        MPI_Isend(&send,Ncities,MPI_INTEGER,destinatario[1],0,MPI_COMM_WORLD,&req1);


      }
      else if(rank==destinatario[1])
      {

        MPI_Recv(&recive,Ncities,MPI_INTEGER,destinatario[0],0,MPI_COMM_WORLD, &stat1);
        p.copy_list(recive,0);//copia la lista delle citta nel primo elemento di popolazioni
        p.order(); // ordina
        MPI_Isend(&send,Ncities,MPI_INTEGER,destinatario[0],0,MPI_COMM_WORLD,&req1);
      }

      if(rank==destinatario[0])
      {
        MPI_Recv(&recive,Ncities,MPI_INTEGER,destinatario[1],0,MPI_COMM_WORLD, &stat1);
        p.copy_list(recive,0);
        p.order();
      }

      if(rank==destinatario[2])
      {

        MPI_Isend(&send,Ncities,MPI_INTEGER,destinatario[3],0,MPI_COMM_WORLD,&req1);


      }
      else if(rank==destinatario[3])
      {

        MPI_Recv(&recive,Ncities,MPI_INTEGER,destinatario[2],0,MPI_COMM_WORLD, &stat1);
        p.copy_list(recive,0);
        p.order();
        MPI_Isend(&send,Ncities,MPI_INTEGER,destinatario[2],0,MPI_COMM_WORLD,&req1);
      }
      if(rank==destinatario[2])
      {
        MPI_Recv(&recive,Ncities,MPI_INTEGER,destinatario[3],0,MPI_COMM_WORLD, &stat1);
        p.copy_list(recive,0);
        p.order();
      }
    }


}while(k<1000);


  for (int i=0; i<Ncities;i++)
  {
    fout1<<posizioni[p.Get(i,0)-1][0]<<" "<<posizioni[p.Get(i,0)-1][1]<<endl;

  }

  fout.close();
  fout1.close();
  fout2.close();
  fout3.close();



  MPI_Finalize();
  return 0;
}

void Init(int rank)
{

	int p1, p2;

  ifstream Primes("Primes");

   if (Primes.is_open())
   {

     int appo1,appo2;
     for(int i=1; i<rank*2;i++)
        Primes >> appo1 >> appo2 ;
      Primes >> p1 >> p2 ;
      cout<<"rank: "<<rank<<"primes  "<<p1<<"  "<<p2<<endl;

	   }
  else cerr << "PROBLEM: Unable to open Primes" << endl;

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
