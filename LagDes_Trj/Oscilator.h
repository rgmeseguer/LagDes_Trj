#pragma once
#include <vector>					//vector
#include <iostream>					//std::cout, std::fixed
#include <math.h>					//sqrt fabs 	
#include <stdlib.h>					//abs srand rand
#include <time.h>					//time	
#include "utilities.h"				//RandomFloat
using namespace std;

class Oscillator
{
	utilities ut;
	typedef double(*PotFunct)(vector<double>, vector<double>);
	PotFunct Pf;
	typedef vector<double>(*GradFunct)(vector<double>, vector<double>);
	GradFunct Gf;

public:

	//creator of the object with the invarian values the coeff and the masses
	Oscillator(vector<double> coeff, vector<double> mass, double(V)(vector<double>, vector<double>), vector<double>(G)(vector<double>, vector<double>));

	vector<double> mu;//Reduced Mass of the 2 oscillators
	size_t size; // Size of the system
	vector<double> c; //coefficients of the system

	vector<double> r; //Position
	vector<double> v; //velocity
	vector<double> ek;// Kinetic Energy
	double ep; // Potential Energy
	vector<double> grad;//Gradient
	vector<double> a;//Acceleration


	//Set the position, the velocities and the Energy of the system
	void setInitP(vector<double>, bool);
	void setInitV(vector<double>, bool);

	//Calculate the properies of the oscillator
	double Etot();// Total energy
	void gradient();//recalculate the gradient
	void acc();//recalculate the Acceleration
	vector<double> p(); //calculate the momenta


	//Finds the point for the bath knowing the position of the remaining DF and the KE being 0
	void keepEnergy0(double Energy);
	//Sets the velocity of the system to keep the energy constant
	bool keepEnergy1(double Energy, int Syst);
	//Sets the position of the Bath to keep the energy Bath has an initial vel
	bool keepEnergy2(double Energy, int Sys, double Psys, double kinbath);
	//Sets the velocity of the two motions to keep the energy constant
	void keepEnergy3(double Energy, int Syst, int Bath);
	//Sets the velocity of the two motions to keep the energy constant I is the index for the velocity of the bath
	void keepEnergy4(double Energy, int Syst, int Bath, int I);



};

