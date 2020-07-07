#include "classi.h"
#include <algorithm>

using namespace std;

population::population ( int Npop_m , int  N_cities_m, Random mrnd , double **posizioni_m)
{
  Npop = Npop_m; //numero di individui (righe)
  N_cities = N_cities_m;// n colonne
  rnd=mrnd;
  posizioni=posizioni_m;
 for(int i = 0; i<Npop; i++)
  {
    individuo ind (N_cities,rnd, posizioni);
    pop.push_back(ind);
    new_pop.push_back(ind);
  }
}

void population::order ()
{
    std::sort(pop.begin(),pop.end());
}


void population::start(){
  for (int i=0; i<N_cities;i++)
  {
    pop[0].Set(i,i+1);

  }

  for (int j=1;j<Npop;j++)
  {

    for (int i=0; i<N_cities;i++)
    {
      pop[j].Set(i,pop[j-1].Get(i));

    }
    int a = floor(rnd.Rannyu(2.,(double)N_cities));
    int b =  floor(rnd.Rannyu(2.,(double) N_cities));


    pop[j].Swap(a,b);

  }
  order();
}
void population::totalcheck()
{
  int sum = 0;
  for (int i=0; i<Npop; i++)
        sum +=pop[i].check( );
  cout<<"check completed:  "<<endl;
  if (sum==Npop)
    cout<<"every individual fulfils the bonds"<<endl;
  else
    cout<<"error with starting sequence"<<endl;
}
double population::Pbc(double r)
{
    r--;
    double box =N_cities - 1; //perchè il primo elemento rimane fermo
    return r - box* floor(r/box)+1;
}

int population::selection()
{
  return (int) (Npop*pow(rnd.Rannyu(),p));
}



bool population::converge()
{
  int sum=0;

  for (int i=0; i<4; i++)
  {

    if(fabs(pop[i].L2()-pop[i+1].L2()) < 0.09)
    {
      sum++;

    }

  }

  if(sum==4)
    return 0;
  else
    return 1;
}

void population::genetic()
{
  order();
  for (int i=0; i<Npop; i+=2)
  {
    individuo i1(N_cities,rnd, posizioni);
    individuo i2(N_cities,rnd, posizioni);
      int a,b;
      do
      {
        a=selection();
        b=selection();

      }while(a ==b);


      i1.copy(pop[a]);
      i2.copy(pop[b]);

      i1.mutation();
      i2.mutation();
      i1.crossover(i2);

      new_pop[i].copy(i1);
      new_pop[i+1].copy(i2);

  }
  for(int i=0;i<Npop;i++)
    pop[i].copy(new_pop[i]);
  order();
}

double population::average()
{
  double ave=0;
  for (int i=0; i<=Npop/2; i++)
    ave += pop[i].L2();
  return ave/((int)Npop/2);

}

//#####################################################################################################################
//#####################################################################################################################

individuo::individuo(int N_cities_m, Random mrnd, double **posizioni_m){

  N_cities=N_cities_m;
  rnd=mrnd;
  posizioni = posizioni_m;

  for (int i=0; i<N_cities; i++)
    ind.push_back(0);
}

bool individuo::check( )
{
  for (int i=0; i<N_cities-1; i++)
  {
    for (int j=i+1;j<N_cities; j++)
    {
      if (ind[i]==ind[j])
      {

        return  0 ;
      }

    }
  }

  return 1;
}

void individuo::mutation( )
{
  int i = floor(rnd.Rannyu(1,4));
  double c = floor(rnd.Rannyu());

  if(c<prob_mutation)
  {


    if(i==1)
    {
      int a = (int ) rnd.Rannyu(1,N_cities-1);
      int b = (int ) rnd.Rannyu(1,N_cities-1);
      std::swap(ind[a],ind[b]);


    }
    if (i==2)
    {

      if ( m2 > N_cities-1 )
        cout<<"ERROR WITH m2 VARIABLE"<< endl;

      int a = floor( rnd.Rannyu(1,N_cities-1));


      for (int j = a; j<=m2+a; j++)
      {

        std::swap(ind[Pbc(j)],ind[Pbc(j+n)]);

      }



    }
    if (i==3)
    {

      if ( m3 > N_cities/2 )
        cout<<"ERROR WITH m3 VARIABLE"<< endl;

      int a =  rnd.Rannyu(1,floor(((N_cities-1)/2)+1));
      int inizio = a;
      int fine = a + m3;

      while ( fine!=inizio || fine != inizio++ )
      {

        std::swap(ind[inizio],ind[fine]);
        inizio = inizio+1;
        fine = fine-1;

      }
    }


  }


}
double individuo::Pbc(double r)
{
    r--;
    double box =N_cities - 1; //perchè il primo elemento rimane fermo
    return r - box* floor(r/box)+1;
}

double individuo::L2 ( )
{


  double dist = 0;

  for (int i=0; i<N_cities-1; i++)
  {

    dist = dist +sqrt (pow(posizioni[ind[i]-1][0]-posizioni[ind[i+1]-1][0],2.)+pow(posizioni[ind[i]-1][1]-posizioni[ind[i+1]-1][1],2.));

  }

  return dist+sqrt(pow(posizioni[ind[N_cities-1]-1][0]-posizioni[ind[0]-1][0],2.)+pow(posizioni[ind[N_cities-1]-1][1]-posizioni[ind[0]-1][1],2.));
}

void individuo::copy(individuo i2)
{
  for (int j=0; j<N_cities; j++)
      ind[j]=i2.Get(j);
}


void individuo::crossover( individuo i2)
{

  double prob = rnd.Rannyu();

  if(prob<prob_crossover)
  {
    int begin =floor( rnd.Rannyu(1 , N_cities-1));
    int end = N_cities-1;



    int *appo1 = new int [N_cities];
    int *appo2 = new int [N_cities];


    int k=begin;
    for (int j=0; j<N_cities; j++)
    {

      for (int i=begin; i<= end; i++ )
      {

        if(ind[i]==i2.Get(j))
        {

          appo1[k]=i2.Get(j);
          k++;

        }
      }
    }

   k=begin;
    for (int j=0; j<N_cities; j++)
    {
      for (int i=begin; i<= end; i++ )
      {
        if(i2.Get(i)==ind[j])
        {

          appo2[k]=ind[j];
          k++;
        }
      }
    }


    for (int i=begin; i<= end; i++ )
    {
      ind[i]=appo1[i];
      i2.Set(i,appo2[i]);
    }
  }



}
