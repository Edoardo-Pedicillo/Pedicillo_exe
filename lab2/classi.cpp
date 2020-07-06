 #include "classi.h"

using namespace std;



double coseno :: Eval ( double rd ) const {

	double x = rd ;



	double z =  M_PI * cos( M_PI * x / 2 )/2;

	return z;
}



double funzione2 :: Eval ( double rd ) const {



	double x =( 1-sqrt(1-rd));


	double z = M_PI* cos(M_PI * x / 2. )/(2*(-2*(x-1)));



	return z;
}
