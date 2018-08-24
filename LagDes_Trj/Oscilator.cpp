#include "Oscilator.h"


Oscillator::Oscillator(vector<double> coeff, vector<double> mass, double (V)(vector<double>, vector<double>), vector<double>(G)(vector<double>, vector<double>))
{
	c = coeff;
	mu = mass;

	size = mu.size();
	vector<double> zerotmp(size, 0.);
	r = zerotmp;
	v = zerotmp;
	a = zerotmp;
	ek = zerotmp;
	grad = zerotmp;
	ep = 0.;
	Pf = V;
	Gf = G;
}

//Set the position, the velocities and the Energies of the system
void Oscillator::setInitP(vector<double> coord, bool Energy)
{
	r = coord;
	grad = Gf(c, r);//the gradient
	for (size_t i = 0; i < size; i++)//the acceleration
	{
		a[i] = -grad[i] / mu[i];
	}

	if (Energy)
	{
		ep = Pf(c, r);//set the potential energy
	}

}
void Oscillator::setInitV(vector<double> veloc, bool Energy)
{
	v = veloc;
	if (Energy)
	{
		//also set the kinetic energy
		for (size_t i = 0; i < size; i++)
		{
			ek[i] = mu[i] / 2 * (v[i] * v[i]);
		}
	}
}

//Calculate the properies of the oscillator
double Oscillator::Etot()
{
	double K = 0.;
	ep = Pf(c, r);//calculates the potential energy

	for (size_t i = 0; i < size; i++)
	{
		ek[i] = mu[i] / 2 * (v[i] * v[i]); // and the kinetic energy
		K += ek[i];
	}

	return ep + K;
}
void Oscillator::gradient()
{
	//recalculate the gradient
	grad = Gf(c, r);
	
}
void Oscillator::acc()//recalculate the Acceleration
{
	grad = Gf(c, r);//the gradient
	for (size_t i = 0; i < size; i++)//the acceleration
	{
		a[i] = -grad[i] / mu[i];
	}
}
vector<double> Oscillator::p()
{
	vector<double> mom;
	for (size_t i = 0; i < size; i++)
	{
		mom.push_back(v[i] * mu[i]);
	}
	return mom;
}


#if 0
//Sets the position of the Bath to keep the energy
bool Oscillator::keepEnergy0(double Energy, int Sys, double Psys)
{
	//We need to searh for the root, the coordinate in the bath that keeps
	//the energy constant with the given momenta and position of the system
	double r0, r1; //Positions
	double dr;//Step size
	double E0, E1;//Energies

#pragma region Initial Conditions
	//Set the initial Velocity from the Momenta of the system
	setInitV({ Psys / mu[Sys] , 0. }, true);
	//Set the initial 3 points to search for the root
	r0 = 1.5;//Beggining
	r1 = 3.0;//End
	dr = (r1 - r0) / 100.;//Divide the range

	//Energy on the Beggining
	vector<double> init = r;
	init.back() = r0;
	setInitP(init, true);
	E0 = Etot() - Energy; //Energy differece

	r1 = r0;//Keep track of the Beggining Point
	E1 = E0;//Keep track of the Beggining EnergyDiff
#pragma endregion

	//Get close to the root
	while (fabs(E0) > 1.e-7)
	{
		//Do another step
		r1 += dr;
		init.back() = r1;
		setInitP(init, true);
		E1 = Etot() - Energy;
		//Did you miss the root?
		if (/*you cross it*/ ((E1*E0 < 0.)
			|| (fabs(E1) > fabs(E0)))
			&& /*and you are close enough to the root*/ (E0 < 10))
		{
			//You missed the root get back and reduce the step
			r1 -= dr; dr *= 0.5;
			init.back() = r1;
			setInitP(init, true);
			E1 = Etot() - Energy;
		}
		else //Reduce the step so it cannot end in a infinite loop
		{
			dr *= 0.999;
		}
		//If no Root do not calculate
		if (fabs(dr) < 1.e-10) {
			return false;
		}
		//Actualize the back position
		E0 = E1;
		r0 = r1;
	}
	init.back() = r0;
	setInitP(init, true);
	return true;
}
//Sets the velocity of the system to keep the energy constant
bool Oscillator::keepEnergy1(double Energy, int Syst)
{
	bool out = true;
	//Set the other velocities to 0
	vector<double> b(size, 0.);
	setInitV(b, true);

	ek[Syst] = Energy - ep;
	if (ek[Syst] < 0.) { out = false; }
	else
	{
		if (fabs(ek[Syst]) < 1.e-6)
		{
			ek[Syst] = 0.;
		}
		v[Syst] = sqrt(2 * ek[0] / mu[0]);
		setInitV(v, true);
	}
	return out;
}
//Sets the position of the Bath to keep the energy but add after the velocity of the bath based on the temperature for Stochastich 
bool Oscillator::keepEnergy2(double Energy, int Sys, double Psys, double kinbath)
{
	//We need to searh for the root, the coordinate in the bath that keeps
	//the energy constant with the given momenta and position of the system
	double r0, r1; //Positions
	double dr;//Step size
	double E0, E1;//Energies

#pragma region Initial Conditions
	//Set the initial Velocity from the Momenta of the system
	setInitV({ Psys / mu[Sys] , 0. }, true);
	//Set the initial 3 points to search for the root
	r0 = 1.5;//Beggining
	r1 = 3.0;//End
	dr = (r1 - r0) / 100.;//Divide the range

						 //Energy on the Beggining
	vector<double> init = r;
	init.back() = r0;
	setInitP(init, true);
	E0 = Etot() - Energy; //Energy differece

	r1 = r0;//Keep track of the Beggining Point
	E1 = E0;//Keep track of the Beggining EnergyDiff
#pragma endregion

	//Get close to the root
	while (fabs(E0) > 1.e-7)
	{
		//Do another step
		r1 += dr;
		init.back() = r1;
		setInitP(init, true);
		E1 = Etot() - Energy;
		//Did you miss the root?
		if (/*you cross it*/ ((E1*E0 < 0.)
			|| (fabs(E1) > fabs(E0)))
			&& /*and you are close enough to the root*/ (E0 < 10))
		{
			//You missed the root get back and reduce the step
			r1 -= dr; dr *= 0.5;
			init.back() = r1;
			setInitP(init, true);
			E1 = Etot() - Energy;
		}
		else //Reduce the step so it cannot end in a infinite loop
		{
			dr *= 0.999;
		}
		//If no Root do not calculate
		if (fabs(dr) < 1.e-10) {
			return false;
		}
		//Actualize the back position
		E0 = E1;
		r0 = r1;
	}
	init.back() = r0;
	setInitP(init, true);//Set the final position of the baht

#pragma region Initial velocity of the bath
	srand(time(NULL));
	int one[2] = { -1,1 };
	int ranIndex = rand() % 2;
	double v01 = one[ranIndex] * sqrt(2 * kinbath / mu[0]);//Velocity of the bath
	setInitV({ v[0],v01 }, true);

#pragma endregion
	return true;
}

#endif // old Keep Energy

//Sets a random velocity for the two motions keeping the energy constant
void Oscillator::keepEnergy3(double Energy, int Syst, int Bath)
{
	srand(time(NULL));
	int one[2] = { -1,1 };
	double kineticTot = Energy - ep;//Calculate the total kinetic energy
	//randomly divide it between the system and the bath
	double kinBath = ut.RandomFloat(0., kineticTot);
	double kinSys = kineticTot - kinBath;

	//cout << kineticTot << " " << kinBath + kinSys << " " << " " << kinBath << " " << kinSys << endl;
	int ranIndex = rand() % 2;
	v[Syst] = sqrt(2 * kinSys / mu[Syst]) * one[ranIndex];

	ranIndex = rand() % 2;
	v[Bath] = sqrt(2 * kinBath / mu[Bath])* one[ranIndex];
	setInitV(v, true);

	cout << Energy << " " << Etot() << " " << ek[0] << " " << ek[1] << endl;

}
//Sets the velocity of the two motions to keep the energy constant I=[0:3999] is the index for the velocity of the bath
void Oscillator::keepEnergy4(double Energy, int Syst, int Bath, int I)
{
#pragma region Select the direction of the system and the bath
	int sysdir, bathdir;
	if (I >= 40000)
	{
		printf("Error:Index too high, Max number of calculations 4000\n");
		terminate();
	}
	else if (I >= 30000)
	{
		I = I - 30000;
		sysdir = -1;
		bathdir = -1;
	}
	else if (I >= 20000)
	{
		I = I - 20000;
		sysdir = 1;
		bathdir = -1;

	}
	else if (I >= 10000)
	{
		I = I - 10000;
		sysdir = -1;
		bathdir = 1;
	}
	else
	{
		sysdir = 1;
		bathdir = 1;
	}

#pragma endregion

	double kineticTot = Energy - ep;//Calculate the total kinetic energy
	double dek = kineticTot / 10000; //Divide the energy into 1000 steps
	double kinSys = dek * I; //Set the KE of the system depending on the index
	double kinBath = kineticTot - kinSys; //Set the KE of the bath

	//cout << kinBath << endl;

	//Set the velocities of the system
	v[Syst] = sqrt(2 * kinSys / mu[Syst]) * sysdir;
	v[Bath] = sqrt(2 * kinBath / mu[Bath])* bathdir;
	setInitV(v, true);

	//cout << Energy << " " << Etot() << " " << ek[0] << " " << ek[1] << endl;

}

//Finds the point for the bath knowing the position of the remaining DF and the KE being 0
void Oscillator::keepEnergy0(double Energy)
{
	//We need to searh for the root, the coordinate in the bath that keeps
	//the energy constant with the given  and position of the system
	double r0, r1; //Positions
	double dr;//Step size
	double E0, E1;//Energies


#pragma region Initial Conditions
				  //Set the initial Velocities to 0
	for (size_t i = 0; i < size; i++) { v[i] = 0.; }
	setInitV(v, true);
	//Set the initial 3 points to search for the root
	r0 = 3.0;//Beggining
	r1 = 1.5;//End
	dr = (r1 - r0) / 100.;//Divide the range

						  //Energy on the Beggining
	vector<double> init = r;
	init.back() = r0;
	setInitP(init, true);
	E0 = Etot() - Energy; //Energy differece

	r1 = r0;//Keep track of the Beggining Point
	E1 = E0;//Keep track of the Beggining EnergyDiff
#pragma endregion

			//Get close to the root
	while (fabs(E0) > 1.e-7)
	{
		//Do another step
		r1 += dr;
		init.back() = r1;
		setInitP(init, true);
		E1 = Etot() - Energy;
		//Did you miss the root?
		if (/*you cross it*/ ((E1*E0 < 0.)
			|| (fabs(E1) > fabs(E0)))
			&& /*and you are close enough to the root*/ (E0 < 10))
		{
			//You missed the root get back and reduce the step
			r1 -= dr; dr *= 0.5;
			init.back() = r1;
			setInitP(init, true);
			E1 = Etot() - Energy;
		}
		else //Reduce the step so it cannot end in a infinite loop
		{
			dr *= 0.999;
		}
		//If no Root do not calculate
		if (fabs(dr) < 1.e-10) { cout << "something went wrong " << r1 << " " << fabs(E0) << endl; break; }
		//Actualize the back position
		E0 = E1;
		r0 = r1;
	}
	init.back() = r0;
	setInitP(init, true);//Set the final position of the bath


}



