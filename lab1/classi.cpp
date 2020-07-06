#include "classi.h"

using namespace std;




expDistribution::expDistribution ( double nlambda )

{

	lambda = nlambda;


};


void expDistribution::SetLambda (double m_lambda)
{

	lambda = m_lambda;
};

double expDistribution::GetLambda ()
{


	return lambda;
};


double expDistribution::Eval ( Random *rd ) const
{




	double y = rd->Rannyu( );



	return -log(1-y)/lambda;

};

LorDistribution::LorDistribution( double ngamma )

{


	gamma = ngamma;

}


void LorDistribution::Setgamma (double m_gamma)
{

	gamma = m_gamma;
};

double LorDistribution::Getgamma ()
{


	return gamma;
};

double LorDistribution::Eval ( Random *rd ) const
{


	double y = rd->Rannyu(  );

	return gamma * tan( 3.14 * (y - 0.5) );


};

Buffon::Buffon ( int sRighe, double sL, double sd, double sxmax, int sM, int sN )
{

	nRighe = sRighe;
	L = sL;
	d = sd;
	xmax = sxmax;
	M = sM;
	N = sN;

};

int Buffon::interseca ( double x1, double y1, double x2, double y2 ) const
{

	double m = (y1 - y2) / (x1 - x2);

	double xriga ;
	double yriga ;

	double distanza ;

	// in teoria questo procedimento non sarebbe propriamente giusto perchè favovirebbe alcune direzione a scapito di altre, ovvero presa una retta e un punto favorisce di più la
	// direzione ( tra le due possibili ) che si trova maggiormente nel quadrato. la cosa si può trascurare se prendiamo un quadrato molto grande.

	if ( x2 < x1) // in questo caso la riga più vicina è quella che approssima per difetto x1
	{
		xriga = trunc ( x1 ) ; // approssima ad intero togliendo la parte decimale, quindi mi permette di calcolare una riga vicina al punto 1 ( le righe sono parallele a y )

		yriga = m * xriga - m * x1 + y1 ; // la y a cui la riga e l'ago si incontrano

		distanza = pow ( (xriga - x1) * ( xriga - x1) + (yriga - y1) * ( yriga - y1 ) , 0.5 ); // distanza x1 punto intersezione



	}
	else
	{
		xriga = trunc ( x1 ) + 1.0;

		yriga = m * xriga - m * x1 + y1;

		distanza = pow ( (xriga - x1) * ( xriga - x1) + (yriga - y1) * ( yriga -y1 ) , 0.5 );



	}

	if ( distanza < L)
		return 1;
	else
		return 0;


};

double Buffon::Eval ( Random *rd ) const {

	double distanza;
	double x1,x2,y1,y2;
	do {
		x1 = rd->Rannyu( 0 , xmax);

		x2 = rd->Rannyu( 0 , xmax);

		y1 = rd->Rannyu( 0 , xmax);

		y2 = rd->Rannyu( 0 , xmax);

	 distanza = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);

	}while (distanza<=L);

	


	return d*this->interseca ( x1, y1, x2, y2) / ( 2 * L ); // Moltiplico per M/N inoltre passo questa espressione che calcola 1/pi in modo tale da
								// poter usare la funzione Simulazione.

}
