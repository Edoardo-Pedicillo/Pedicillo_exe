#define _USE_MATH_DEFINES //per usare pi greco
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

using namespace std;


class FunzioneBase {

	public:

	virtual double Eval (double r []) const =0;

};


class psi_1s : public FunzioneBase
{

	public:

	psi_1s ( ) {;} ;

	double Eval (double r []) const;



};

class psi_2p : public FunzioneBase
{

	public:

	psi_2p ( ) {;} ;

	double Eval (double r []) const;



};
