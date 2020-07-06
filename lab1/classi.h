#define _USE_MATH_DEFINES //per usare pi greco
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

using namespace std;

// Definisco la classe astratta FunzioneBase e le relative classi derivate poichè l'es 1.1.1 e 1.1.2 hanno lo stesso codice ad eccezione del termine che si somma nella variabile sum
//in questo modo evito di ripetere parti di codice

class FunzioneBase {

	public:

	virtual double Eval (double r) const =0; // virtual: il metodo può essere sovrascritto dalle classi figle
						//  = 0 vuol dire che questo metodo non può essere implementato nella classe stessa
						
	virtual double Eval (Random *rd) const = 0;

};


class unifDistribution : public FunzioneBase // Semplicemente restituisce il numero Random
{

	public:

	unifDistribution ( ) {;} ;

	double Eval (double r ) const { return r; };

	double Eval (Random *rd) const { return rd->Rannyu( ) ; };


};


class dev : public FunzioneBase {


	public:

	double Eval (double r ) const {return (r-0.5)*(r-0.5);};
	double Eval (Random *rd) const { return 0;};


};

class expDistribution : public FunzioneBase {

	public:

	expDistribution ( double nlambda);


	void SetLambda (double m_lambda);
	double GetLambda () ;
	double Eval ( Random *rd ) const ;
	double Eval (double r ) const {return 0;};

	private:

	double lambda;



};


class LorDistribution: public FunzioneBase {

	public:

	LorDistribution( double ngamma);


	void Setgamma (double m_gamma);
	double Getgamma () ;
	double Eval ( Random *rd ) const ;
	double Eval (double r ) const {return 0;};

	private:

	double gamma;




};

class Buffon : public FunzioneBase {

	public:

	Buffon ( int sRighe, double sL, double sd, double sxmax , int sM, int sN );

	double Eval ( Random *rd ) const ;
	double Eval (double r ) const {return 0;};

	int interseca ( double x1, double y1, double x2, double y2 ) const ;


	private:

	int nRighe;
	double L;
	double d;
	double xmax;
	int M;
	int N;
};
