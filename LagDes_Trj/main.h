#pragma once
//#include <cstdio>
#include <vector>				// std::vector
#include <iostream>				// std::cout, std::fixed
#include <iomanip>				// std::setprecision
#include <fstream>				// std::ofstream
#include <stdlib.h>				// srand, rand 
#include <time.h>				// time	
#include <sstream>				// stringstream

#include "utilities.h"			// Utilities fucntion
#include "Oscilator.h"			// Oscilator class
#include "LagDesc.h"			// LD calculator
#include "Dynamic.h"			// Dynamics calculator
#include "Point.h"				// Surface calculator

using namespace std;
utilities ut;


//Potential and Gradient of the system
double DOS_V(vector<double> coeff, vector<double> r)
{
	
	double ep = 0.;
	for (int i = 0; i < 5; ++i)
	{
		ep += coeff[i] * ut.powerd(r[0], i);
	}
	ep += coeff[5] * ut.powerd(r[1] - coeff[6], 2);
	ep += coeff[7] / ut.powerd(r[1] - r[0], 12);

	return ep;
}
vector<double> DOS_G(vector<double> coeff, vector<double> r)
{
	
	vector<double> grad(2, 0.);
	grad[0] = 0;
	for (int j = 1; j < 5; ++j)
	{
		grad[0] += double(j)*coeff[j] * ut.powerd(r[0], j - 1);
	}
	grad[0] += 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	grad[1] = 0;
	grad[1] = 2. * coeff[5] * (r[1] - coeff[6]) - 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	return grad;
}

//dot momenta for the LD calculation
vector<double> pdot(LagDesc LD, Oscillator osc, vector<double> R0)
{
	vector<double> p = { 0.,0. };
	p[1] = 2 * R0[1] * osc.c[5] - 2 * osc.c[5] * osc.r[2];

	for (size_t i = 0; i < 4; i++)
	{
		p[0] += LD.D[i] * ut.powerd(osc.r[0], i);
	}

	return p;
}
