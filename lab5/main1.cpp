#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.h"


using namespace std;

int main () {

FunzioneBase *f = new psi_1s (  );
FunzioneBase *g = new psi_2p ( );
Random *rnd = new Random ();
inizRandom ( rnd );
int acc=0;
double x[3];
double delta =1;

x[0]=0.;
x[1]=0.;
x[2]=0.;

int Nthrow = 1000000;
int Nblocks = 100;
double ave [Nblocks];
double av2 [Nblocks];
double sum_prog [Nblocks];
double su2_prog [Nblocks];
double err_prog [Nblocks];
double rs[Nthrow];
// I put the starting point far from the origin

x[0]=100.;
x[1]=100.;
x[2]=100.;

ofstream fout ("Positions100.dat");

for (int i=0; i<1000; i++)
{
	acc += Metropolis (x,f,delta,rnd,0);
	fout<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
}


fout.close();

// Starting point in (0,0,0)

x[0]=0;
x[1]=0;
x[2]=0;

ofstream fout0 ("Positions0.dat");

for (int i=0; i<5000; i++)
{
	acc += Metropolis (x,f,delta,rnd,0);

	fout0<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
}


fout0.close();

//:::::::::::::::::::::UNIFORM:::::::::::::::::::
cout<<":::::Uniform Distribution:::::"<<endl;

// psi_1s
delta=1.2;
x[0]=0;
x[1]=0;
x[2]=0;

//Equilibration

for (int i=0; i<100; i++)
{
		acc += Metropolis (x,f,delta,rnd,0);
}

acc=0;

// <r>

for (int i=0; i<Nthrow; i++)
{
	acc += Metropolis (x,f,delta,rnd,0);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	rs[i] = r;
}
cout<<"psi_1s:  "<<endl;
cout<<"acceptance:  "<<(double)acc/Nthrow<<endl;

Simulazione ( Nthrow,Nblocks, rs, ave, av2, sum_prog, su2_prog, err_prog );

ofstream fout1 ("psi1s_unif.dat");

for (int i=1; i<Nblocks+1; i++)
{


	fout1 << i <<" "<<sum_prog[i-1] <<" "<<err_prog[i-1] << endl;



}
fout1.close();

//psi_2p

delta = 3;
acc = 0;
ofstream foutb ("Positions0-2p.dat");
x[0]=1;
x[1]=1;
x[2]=1;

for (int i=0; i<100; i++)
{
		acc += Metropolis (x,g,delta,rnd,0);
}

acc=0;

for (int i=0; i<Nthrow; i++)
{
	acc += Metropolis (x,g,delta,rnd,0);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	rs[i] = r;
	if(i<5000)
		foutb<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
}

cout<<"psi_2p:   "<<endl;
cout<<"acceptance:  "<<(double)acc/Nthrow<<endl;

Simulazione ( Nthrow,Nblocks, rs, ave, av2, sum_prog, su2_prog, err_prog );

ofstream fout2 ("psi2p_unif.dat");

for (int i=1; i<Nblocks+1; i++)
{


	fout2 << i <<" "<<sum_prog[i-1] <<" "<<err_prog[i-1] << endl;



}
fout2.close();

//:::::::::::::::::::::Gaussian:::::::::::::::
cout<<endl<<":::::Gaussian Distribution:::::"<<endl;

//psi 1s

x[0]=0.;
x[1]=0.;
x[2]=0.;
delta= 0.75;
for (int i=0; i<100; i++)
{
		acc += Metropolis (x,f,delta,rnd,1);
}

acc=0;

ofstream fout00 ("Positions0G.dat");

for (int i=0; i<5000; i++)
{
	acc += Metropolis (x,f,delta,rnd,1);

	fout00<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
}


fout00.close();

acc=0;

for (int i=0; i<Nthrow; i++)
{
	acc += Metropolis (x,f,delta,rnd,1);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	rs[i] = r;
}
cout<<"psi_1s:  "<<endl;
cout<<"acceptance:  "<<(double)acc/Nthrow<<endl;

Simulazione ( Nthrow,Nblocks, rs, ave, av2, sum_prog, su2_prog, err_prog );

ofstream fout3 ("psi1s_gauss.dat");

for (int i=1; i<Nblocks+1; i++)
{


	fout3 << i <<" "<<sum_prog[i-1] <<" "<<err_prog[i-1] << endl;



}
fout3.close();

//psi_2p

x[0]=1.;
x[1]=1.;
x[2]=1.;
delta = 1.9;

for (int i=0; i<100; i++)
{
		acc += Metropolis (x,f,delta,rnd,1);
}

acc = 0;

for (int i=0; i<Nthrow; i++)
{
	acc += Metropolis (x,g,delta,rnd,1);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	rs[i] = r;
}

cout<<"psi_2p:   "<<endl;
cout<<"acceptance:  "<<(double)acc/Nthrow<<endl;

Simulazione ( Nthrow,Nblocks, rs, ave, av2, sum_prog, su2_prog, err_prog );

ofstream fout4 ("psi2p_gauss.dat");

for (int i=1; i<Nblocks+1; i++)
{


	fout4 << i <<" "<<sum_prog[i-1] <<" "<<err_prog[i-1] << endl;



}
fout4.close();


	return 0;
}
