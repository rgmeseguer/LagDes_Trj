#include "LagDesc.h"
LagDesc::LagDesc()
{
}


LagDesc::~LagDesc()
{
}

#if KEY_M1
double LagDesc::M1(vector<double> v)
{
	double sqV = 0.;
	for (size_t i = 0; i < v.size(); i++)
	{
		sqV += (v[i] * v[i]);
	}
	return sqrt(sqV);
}
#endif
#if KEY_M5
double LagDesc::M5(vector<double> v, vector<double> a)
{
	double sqV = 0.;
	double sqA = 0.;
	double sqB = 0.;
	for (size_t i = 0; i < v.size(); i++)
	{
		sqV += (v[i] * v[i]);
		sqA += (a[i] * a[i]);
		sqB += (a[i] * v[i]);
	}
	
	return 1. / (1. + sqrt(ut.absol((sqV*sqA - sqB*sqB) / ut.powerd(sqV, 3))));
}
#endif

#if KEY_MS
void LagDesc::LoadMs(Oscillator osc)
{
	D[0] = osc.r[0] * 2 * osc.c[2];
	D[1] = 6 * osc.r[0] * osc.c[3] - 2 * osc.c[2];
	D[2] = 12 * osc.r[0] * osc.c[4] - 6 * osc.c[3];
	D[3] = -12 * osc.c[4];

}

vector<double> LagDesc::Ms(vector<double> rdot, vector<double> pdot)
{
	vector<double> M(rdot.size());

	for (size_t i = 0; i < rdot.size(); i++)
	{
		M[i] = sqrt(pdot[i] * pdot[i]) + sqrt(rdot[i] * rdot[i]);
	}
	return M;
}
#endif

#if KEY_MF

vector<double> LagDesc::Mf(vector<double> momenta, vector<double> rdot)
{
	vector<double> M(rdot.size());

	for (size_t i = 0; i < rdot.size(); i++)
	{
		M[i] = sqrt( momenta[i] * rdot[i] );
	}
	return M;
}
#endif
