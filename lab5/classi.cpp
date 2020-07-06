#include "classi.h"

using namespace std;





double psi_1s::Eval ( double x[] ) const
{
	double r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	return exp(-2*r)/M_PI;

};

double psi_2p::Eval ( double x[] ) const
{
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double cos_teta = x[2]/r;
	return exp(-r)* (2/M_PI)*r*r*cos_teta*cos_teta/64 ;

};
