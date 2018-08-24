#pragma once
#include <stdlib.h>				//size_t
#include <math.h>				//sqrt
#include <vector>				//vector	
using namespace std;

//Constant definitions
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)


class utilities
{
public:
	utilities();

	// Extended version of power function that can work
	//for double x and negative y
	double powerd(double x, int y);

	//Random float betwen 2 Floats
	float RandomFloat(float a, float b);

	//Creates a Matrix of size = size1 X size2
	double** create_matrix(size_t size1, size_t size2);
	
	//Box Muller Random Numbers
	//The Box-Muller transform takes two random variables,
	//evenly distributed in the interval (0,1) and transforms them to two independent deviates,
	//which are sampled from a Gaussian distribution. 
	float ran2(long *idum);
	vector<double> box_muller(int nstep, long *ran);
	//Returns the absolute value of a double
	double absol(double val);
};

