#pragma once
#include <iostream>					//cout	
#include <sstream>					//stringstream
#include <iomanip>					//setprecision
#include <fstream>					//ofstream
//#include <stdio.h>

#include "Oscilator.h"

using namespace std;

#pragma once
class Dynamic
{
	double T=0;
public:
	int nstep = 0;
	double dt = 1;
	double dt2 = dt/2.;
	//Langevin Variables
	double A, B;
	double sigma;
	vector<double> gasdev;
	double kinbath;

	//Constructor of the dynamic
	Dynamic();
	//Set the parameters of the Dynamic
	void setTime(double t);
	void setTimeStep(double t);
	//Set the Langevin properties
	void setLangevin(double mu, double beta, double gamma);
	void saveTrj(Oscillator &osc, ofstream &sfile);

	
#if KEY_DETERM
	//Deterministic Relax Dynamic TimeStep
	void DStepRelax(Oscillator & osc);
	//Deterministic Production Dynamic TimeStep
	void DynamicStep(Oscillator&);
	//Calculate a Full Deterministic Trj
	void DynamicTrj(Oscillator& osc, bool print, bool relax);
	vector<vector<double>> ProdTrj(Oscillator & osc, int T);

#else
	//Langevin Dynamic Relax Time Step
	void DStepRelax(Oscillator & osc, /*double B, double A,*/ double rtherm);
	//Langevin Dynamic Production Time Step
	void DynamicStep(Oscillator & osc, /*double B, double A,*/ double rtherm);
	//Calculate a Full Langevin Trj
	void DynamicTrj(Oscillator& osc, bool print, bool relax);

#endif // DETERM








};

