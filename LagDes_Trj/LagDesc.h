#pragma once
#include <vector>					//vector	
#include <math.h>					//sqrt

#include "utilities.h"				//powerd
#include "Oscilator.h"				//Oscillator
using namespace std;

class LagDesc
{
	utilities ut;
public:
	LagDesc();
	~LagDesc();
	
	vector<double> D = { 0,0,0,0 };

	double M1(vector<double> p);
	double M5(vector<double> v, vector<double> a);

	void LoadMs(Oscillator osc);
	vector<double> Ms(vector<double> rdot, vector<double> pdot);
	vector<double> Mf(vector<double> momenta, vector<double> rdot);


};

