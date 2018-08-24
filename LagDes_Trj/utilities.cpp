#include "utilities.h"

utilities::utilities()
{
}

//Random float betwen 2 Floats
float utilities::RandomFloat(float a, float b)
{
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

// Extended version of power function that can work
//for double x and negative y
double utilities::powerd(double x, int y)
{
	double temp;
	if (y == 0)
		return 1;
	temp = powerd(x, y / 2);
	if ((y % 2) == 0) {
		return temp * temp;
	}
	else {
		if (y > 0)
			return x * temp * temp;
		else
			return (temp * temp) / x;
	}
}

//Creates a Matrix of size = size1 X size2
double** utilities::create_matrix(size_t size1, size_t size2)
{
	double** m = new double*[size1];
	for (size_t i = 0; i < size1; ++i)
	{
		m[i] = new double[size2];
	}
	return m;
}

//Box Muller Random Numbers
//The Box-Muller transform takes two random variables,
//evenly distributed in the interval (0,1) and transforms them to two independent deviates,
//which are sampled from a Gaussian distribution. 
float utilities::ran2(long *idum)
/*Long period(> 2 × 10 18) random number generator of L’Ecuyer with Bays - Durham shuffle
and added safeguards.Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values).Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence.RNMX should approximate the largest floating
value that is less than 1.*/
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0) {								//Initialize.
		if (-(*idum) < 1) *idum = 1;				//Be sure to prevent idum = 0.
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) {			//Load the shuffle table(after 8 warm - ups).
			k = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;								//Start here when not initializing.
	*idum = IA1 * (*idum - k * IQ1) - k * IR1;			//Compute idum = (IA1*idum) % IM1 without
	if (*idum < 0) *idum += IM1;					//	overflows by Schrage’s method.
	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;			//Compute idum2 = (IA2*idum) % IM2 likewise.
	if (idum2 < 0) idum2 += IM2;
	j = iy / NDIV;									//Will be in the range 0..NTAB - 1.
	iy = iv[j] - idum2;								//Here idum is shuffled, idum and idum2 are
	iv[j] = *idum;									//	combined to generate output.
	if (iy < 1) iy += IMM1;							//Because users don’t expect endpoint values.
	if ((temp = AM * iy) > RNMX) return RNMX;
	else return temp;
}
vector<double> utilities::box_muller(int nstep, long *ran)
{
	double u1, u2;
	double term1;
	vector<double> gasdev(nstep, 0.);
	for (size_t i = 0; i < nstep; i += 2)
	{
		u1 = ((double)rand() / (RAND_MAX));
		u2 = ((double)rand() / (RAND_MAX));

		u1 = ran2(ran);
		u2 = ran2(ran + 1);

		term1 = sqrt(-2 * log(u1));
		gasdev[i] = term1 * cos(2 * M_PIl * u2);
		gasdev[i + 1] = term1 * sin(2 * M_PIl * u2);
	}
	return gasdev;
}

//Returns the absolute value of a double
double utilities::absol(double val)
{
	if (val < 0.) { return val * -1; }
	else{ return val; }
}
