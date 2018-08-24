#include "main.h"

int main(int argc, char *argv[])
{
	if ((argc!=5))													//Check the correct Parameters
	{
		cout << "Program Usage" << endl;
		cout << "LagD_Trj.out [tau] [bathMass] [TrjTime > (2xtau+1)] [Trj Number]" << endl;
		return 1;
	}
	srand(time(NULL));												//Initialize the random seed


#pragma region Input variables
	double tau = strtof(argv[1], NULL);								//Tau value
	double bathMass = strtof(argv[2], NULL);						//Mass of the Bath
	double TrjTime = strtof(argv[3], NULL);							//Time of the Trj
	int I = int(strtof(argv[4], NULL));								//Number of the trajectory
#pragma endregion
#pragma region Variables for printing
	string BthM, stau;												//String version of tau and Mass
	stringstream stream;
	ofstream outputTrj;												//File of the trj

	//Print the bath mass in the correct format 
	if (bathMass < 1) { stream << fixed << setprecision(1) << bathMass;	BthM = stream.str(); }
	else { stream << fixed << setprecision(0) << bathMass; BthM = stream.str(); }
	stream.str(string());											//clear the stringstream
	//Print the Tau in the correct format
	if (tau < 0) { stream << fixed << setprecision(1) << tau; stau = stream.str(); }
	else { stream << fixed << setprecision(0) << tau; stau = stream.str(); }
	stream.str(string());											//clear the stringstream

#if KEY_DETERM==1
	string calc = "LDT" + stau + "m" + BthM;							//Name of the calculation
#else
	string calc = "sLDT" + stau + "m" + BthM;							//Name of the calculation
#endif

#if KEY_SAVETRJ==1
	outputTrj.open(calc + "_" + to_string(I) + ".trj", ios::out | ios::trunc); //File of the trj
#endif

#pragma endregion
#pragma region System Initial Conditions and Variables
	//Set Initial Conditions
	vector<double> coeff = { 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 };	//Coefficients from the oscillator
	vector<double> mu = { 1.,bathMass };																//Mass Values
	vector<double> R0 = { 1.36561 ,2.161769 };															//Initial Position
	double dt = 1.e-3;																					//time precission
	double Energy = 3.691966889;																		//Energy of the system

	Point Surface;													//Create the surface of the where the LD values will be saved
	Oscillator Osc(coeff, mu, DOS_V, DOS_G);						//Initiate the Oscillator
	Osc.setInitP(R0, true);											//Initiate at the shadle point
	//Osc.keepEnergy3(Energy, 0, 1);									//Set initial Random Velocities to keep energy
	Osc.keepEnergy4(Energy, 0, 1, I);								//Set initial Random Velocities to keep energy
	Dynamic Dyn;													//and the dynamics
	Dyn.setTimeStep(dt);											//the Time step
	Dyn.setTime(TrjTime);											//and the total time (nsteps = totalTime/timeTtep)

#if KEY_DETERM==0													//Variables fot the Stochastic Dynamics
	double beta = 4.;												//Effective temperature
	double gamma = 5.;												//Disipation factor
	Dyn.setLangevin(Osc.mu.back(), beta, gamma);

	long T = time(0);
	vector<double> gasdev = ut.box_muller(Dyn.nstep, &T);			//Gas deviation
	double rtherm;
#endif

#pragma endregion
#pragma region Tool Variables
	LagDesc LD;														//Lagrangian Descriptor variable
	vector<string> key = { "","","","" };							//Variables need for create the Points surface
	string sr0, sr1, sp0, sp1;

	double r0, r1, p0, p1;											//Initial conditions in the surface, position r0 r1 and momentum p0 p1
	vector<vector<string>> OpenPoints;								//Tracker of the points at which the LD is beig calculated							
																	
	vector<double> M1, M5;											//Buffer of the LD value at each step
	vector<vector<double>> Ms, Mf;									//Buffer of the LD value at each step
	vector<double> M = { 0.,0. };									//Buffer of the LD value
	vector<double> tmp = { 0.,0. };

#if KEY_MS==1
	LD.LoadMs(Osc);													//Load the MS LD
#endif

#pragma endregion

	for (size_t j = 0; j < Dyn.nstep; j++)							//Start the trajectory	
	{
#pragma region Dynamic Step
#if KEY_DETERM==1
		Dyn.DynamicStep(Osc);										//Dynamic step
		//cout << Osc.Etot() << endl;
#else
		rtherm = Dyn.sigma*gasdev[j];								//Bath Random
		Dyn.DynamicStep(Osc, rtherm);								//Dynamic step
#endif
#if KEY_SAVETRJ==1
		Dyn.saveTrj(Osc, outputTrj);								//Save the trajectory
#endif

#pragma endregion
#pragma region Lagrangian Descriptor Calculation
																	//
#if KEY_MS==1
		tmp = LD.Ms(Osc.v, pdot(LD, Osc, R0));
		tmp[0] *= dt;
		tmp[1] *= dt;
		Ms.push_back(tmp);											//Save it
		M[0] += Ms.back()[0];										//Add it
		M[1] += Ms.back()[1];
#endif
#if KEY_MF==1
		tmp = LD.Mf(Osc.p(), Osc.v);
		tmp[0] *= dt;
		tmp[1] *= dt;
		Mf.push_back(tmp);											//Save it
		M[0] += Mf.back()[0];										//Add it
		M[1] += Mf.back()[1];
#endif

#pragma endregion

		if (j*Dyn.dt >= tau)
		{
#pragma region Actualize Point
			for (size_t k = 0; k < OpenPoints.size(); k++)			//Actualize every opened point
			{
				key = OpenPoints[k];								//Get the Key of an open point
				Surface.length[key] += dt;							//Increase the length of the point
																	//Add the new LD value 
#if KEY_MS==1
				Surface.P[key].back()[0] += Ms.back()[0];
				Surface.P[key].back()[1] += Ms.back()[1];
#endif
#if KEY_MF==1
				Surface.P[key].back()[0] += Mf.back()[0];
				Surface.P[key].back()[1] += Mf.back()[1];
#endif

			}

#pragma endregion
#pragma region Key Write
			r0 = Osc.r[0];											//Get the values for the new key
			r1 = Osc.r[1];											//Get the values for the new key
			p0 = Osc.v[0] / Osc.mu[0];
			p1 = Osc.v[1] / Osc.mu[1];
			//Print the new key in correct string format
			stream << fixed << setprecision(3) << r0; sr0 = stream.str(); stream.str(string());
			stream << fixed << setprecision(3) << p0; sp0 = stream.str(); stream.str(string());
			stream << fixed << setprecision(3) << r1; sr1 = stream.str(); stream.str(string());
			stream << fixed << setprecision(7) << p1; sp1 = stream.str(); stream.str(string());

			key = { sr0, sp0, sr1, sp1 };							//The location of the new point


#pragma endregion
#pragma region Opening and Closing Points
																	//Open new point only until reach Time-tau 
			if (j*Dyn.dt < Dyn.nstep*Dyn.dt - tau)
			{
				//If the key exists and the point is not open, open the point again				
				if (Surface.isKey(key))
				{
					if (Surface.complete[key] == false)
					{

					}
					else
					{
						Surface.open(key);
						Surface.P[key].push_back(M);
						Surface.length[key] = tau;
						OpenPoints.push_back(key);
					}

				}
				//If it does not exist create a new one		
				else
				{
					Surface.addPoint(key);
					Surface.open(key);
					Surface.P[key].push_back(M);
					Surface.length[key] = tau;
					OpenPoints.push_back(key);
				}
			}
			//If the oldest open point is finished then close it
			if (Surface.length[OpenPoints[0]] >= (2 * tau))
			{
				Surface.close(OpenPoints[0]);
				OpenPoints.erase(OpenPoints.begin());
			}

#pragma endregion
#pragma region Clean Buffer										
			//Clean the first saved point 
#if KEY_MS==1
			M[0] -= Ms[0][0];
			M[1] -= Ms[0][1];
			Ms.erase(Ms.begin());
#endif
#if KEY_MF==1
			M[0] -= Mf[0][0];
			M[1] -= Mf[0][1];
			Mf.erase(Mf.begin());
#endif

#pragma endregion
		}

	}
#pragma region Printing
	outputTrj.close();									//Close The saving TRJ
	outputTrj.open(calc + "_" + to_string(I) + "_" + "LD.txt", ios::out | ios::trunc);
	Surface.SavePav(outputTrj);											//Calc the average and print it
	outputTrj.close();
#pragma endregion

	return 0;
}