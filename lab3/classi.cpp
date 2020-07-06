#include "classi.h"

using namespace std;



PriceCall::PriceCall ( double m_S, double m_T, double m_K, double mr, double msigma, double mtheta)
{
	S_0 = m_S;
	T = m_T;
	K = m_K;
	r = mr;
	sigma = msigma;
	theta = mtheta;
}

double PriceCall:: Eval (Random *rd) const
{
	double Z = rd->Gauss(0,1);
	double S = S_0 * exp( ( r - sigma * sigma / 2 ) * T + sigma * Z * pow( T , 0.5 ) );
	double C = exp ( - r*T ) * max(0. , theta *( S-K) );
	return C;
}

PriceCall2::PriceCall2 ( double m_S, double m_T, double m_K, double mr, double msigma, int nstep, double mtheta )
{
	S_0 = m_S;
	T = m_T;
	K = m_K;
	r = mr;
	sigma = msigma;
	stepT = nstep;
	theta = mtheta;
}

double PriceCall2 :: Eval (Random *rd) const
{
	double S = S_0;
	double C;
	for (int i=1; i<stepT+1; i++)
	{
		double Z = rd->Gauss(0,1);

		S = S * exp( ( r - sigma * sigma / 2 ) * (T/stepT)+ sigma * Z * pow( T/stepT , 0.5 ) );

	}
	C = exp ( - r*T ) * max(0. , theta * (S-K ));
	return C;
}
