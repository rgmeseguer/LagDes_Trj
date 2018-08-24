#include "Dynamic.h"

//Set the parameters of the Dynamic
Dynamic::Dynamic(){}
void Dynamic::setTime(double t)
{
	T = t;
	nstep = abs(int(t / dt));
}
void Dynamic::setTimeStep(double t) { dt = t; dt2 = dt / 2.;}
void Dynamic::saveTrj(Oscillator &osc, ofstream &sfile)
{
	double r;
	double p;
	stringstream stream;
	string sr, sp;

	for (size_t i = 0; i < osc.size; i++)
	{
		r = osc.r[i];
		p = osc.v[i] / osc.mu[i];

		stream << fixed << setprecision(3) << r;				//Print the values in correct string format
		sr = stream.str(); stream.str(string());
		stream << fixed << setprecision(5) << p;				//Print the values in correct string format
		sp = stream.str(); stream.str(string());
		sfile << sr << ' ' << sp << ' ';
	}
	sfile << endl;

}



#if KEY_DETERM
//Deterministic Relax Dynamic
void Dynamic::DStepRelax(Oscillator & osc)
{
	osc.v.back() += osc.a.back() * dt2;
	osc.r.back() += osc.v.back() * dt;
	osc.acc();
	osc.v.back() += osc.a.back() * dt2;
}
//Deterministic Production Dynamic
void Dynamic::DynamicStep(Oscillator & osc)
{
	for (int unsigned i = 0; i < osc.size; i++)
	{
		//osc.v[i] += osc.a[i] * dt2;
		osc.r[i] += (osc.v[i] += osc.a[i] * dt2) * dt;
	}
	osc.acc();
	for (int unsigned i = 0; i < osc.size; i++)
	{
		osc.v[i] += osc.a[i] * dt2;
	}
}
//Calculate a Full Deterministic Trj
void Dynamic::DynamicTrj(Oscillator & osc, bool print, bool relax)
{
	double Eerror, E0 = osc.Etot();

	if (relax)
	{
		for (int i = 0; i < nstep; i++)
		{
			//Perfoms a Dynamic Step and prints it
			DStepRelax(osc);
			Eerror = fabs((osc.Etot() - E0) / osc.Etot());
			//Print the results
			if (print) { cout << i + 1 << ' ' << osc.r[0] << ' ' << osc.r[1] << ' ' << Eerror << endl; }
		}

	}
	else
	{
		for (int i = 0; i < nstep; i++)
		{
			//Perfoms a Dynamic Step and prints it
			DynamicStep(osc);
			Eerror = fabs((osc.Etot() - E0) / osc.Etot());
			//Print the results
			if (print) { cout << i + 1 << ' ' << osc.r[0] << ' ' << osc.r[1] << ' ' << Eerror << endl; }
		}

	}
}
vector<vector<double>> Dynamic::ProdTrj(Oscillator & osc, int T)
{
	//double Eerror, E0 = osc.Etot();
	//long tm = time(0);
	int N = 100;

	int j = 0;
	vector<vector<double>> rout;
	vector<double> zerotmp(2, 0.);
	int Tstep = N * T;
	for (int i = 0; i < Tstep; i++)
	{
		if (j == N)
		{
			zerotmp[0] = osc.r[1];
			zerotmp[1] = osc.v[1];
			rout.push_back(zerotmp);  j = 0;
		}
		//Perfoms a Dynamic Step and prints it
		DStepRelax(osc);
		//Eerror = fabs((osc.Etot() - E0) / osc.Etot());
		//Print the results
		j += 1;
	}
	return rout;
}

#else //Langevin Dynamic
void Dynamic::setLangevin(double mu /*bath mass*/, double beta /*Effective Temperature*/, double gamma /*Dissipation Factor*/)
{
	sigma = sqrt(2 * gamma*fabs(dt) / beta);
	A = gamma * fabs(dt) / (2 * mu);
	B = 1 / (1 + A);
	A = (1 - A) / (1 + A);
	kinbath = 1 / (2 * beta);
}

//Langevin Relax Dynamic
void Dynamic::DStepRelax(Oscillator & osc, /*double B, double A,*/ double rtherm)
{
	//Performs a single molecular dynamics time step via the algorithm of Grønbech-Jensen & Farago
	//initial acceleration
	//vector<double> a0 = osc.a;
	vector<double> g0 = osc.grad;

	//Set the new position
	//Bath has to be always the last defined degree of freedom
	osc.r.back() += B * dt*(osc.v.back() - (osc.grad.back() * dt + rtherm) / (2 * osc.mu.back()));
	osc.setInitP(osc.r, true);
	//Set the new Velocity
	//Bath has to be always the last defined degree of freedom
	osc.v.back() = A * osc.v.back() - dt * (A*g0.back() + osc.grad.back()) / (2 * osc.mu.back()) + B * rtherm / osc.mu.back();
	osc.setInitV(osc.v, true);

}
//Langevin Production Dynamic
void Dynamic::DynamicStep(Oscillator & osc, /*double B, double A,*/ double rtherm)
{
	//Performs a single molecular dynamics time step via the algorithm of Grønbech-Jensen & Farago
	//initial acceleration
	vector<double> a0 = osc.a;
	//Set the new position
	for (int unsigned i = 0; i < osc.size - 1; i++)
	{
		osc.r[i] += dt * (osc.v[i] + osc.a[i] * dt2);
	}
	//Bath has to be always the last defined degree of freedom
	osc.r.back() += B * dt*(osc.v.back() + osc.a.back() * dt2 + rtherm / (2 * osc.mu.back()));
	osc.acc();
	//Set the new Velocity
	for (int unsigned i = 0; i < osc.size - 1; i++)
	{
		osc.v[i] += dt2 * (a0[i] + osc.a[i]);
	}
	//Bath has to be always the last defined degree of freedom
	osc.v.back() = A * osc.v.back() + dt2 * (a0.back() + osc.a.back()) + B * rtherm / osc.mu.back();
}
//Calculate a Full Langevin Trj
void Dynamic::DynamicTrj(Oscillator & osc, bool print, bool relax)
{
	//double Eerror, E0 = osc.Etot();
	double rtherm;

	if (relax)
	{
		for (int i = 0; i < nstep; i++)
		{
			rtherm = sigma * gasdev[i];
			//Perfoms a Dynamic Step and prints it
			DStepRelax(osc,/* B, A,*/ rtherm);
			//Eerror = fabs((osc.Etot() - E0) / osc.Etot());
			//Print the results
			if (print) { cout << i + 1 << ' ' << osc.r[0] << ' ' << osc.r[1] << ' ' << osc.Etot() << endl; }
		}
	}
	else
	{
		for (int i = 0; i < nstep; i++)
		{
			//Perfoms a Dynamic Step and prints it
			rtherm = sigma * gasdev[i];
			DynamicStep(osc, /*B, A,*/ rtherm);
			//Eerror = fabs((osc.Etot() - E0) / osc.Etot());
			//Print the results
			if (print) { cout << i + 1 << ' ' << osc.r[0] << ' ' << osc.r[1] << ' ' << osc.Etot() << endl; }
		}
	}
}

#endif // DETERM




