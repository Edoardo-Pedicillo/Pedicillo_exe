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
						// Eval prende come parametro un double che per gli es 1.1.1 e 1.1.2 è un numero random generato
	virtual double Eval (Random *rd) const = 0;

};



class coseno: public FunzioneBase
{

	public:
		coseno ( ) {;} ;

		double Eval ( Random *rd ) const {return 0;};
		double Eval (double r ) const ;



};


class funzione2 : public FunzioneBase
{

	public:
		funzione2 ( ) { ; } ;

		double Eval ( Random *rd ) const {return 0;};
		double Eval (double r ) const ;


};



class identita : public FunzioneBase
{
	public:

	identita() { ; } ;

	double Eval ( Random *rd ) const {return 0;};
	double Eval (double r ) const {return r;};

};
