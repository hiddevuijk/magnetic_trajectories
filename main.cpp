
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "vecmanip.h"
#include "init_random.h"
#include "deriv_mass.h"
#include "deriv_od.h"
#include "integrate.h"

#include <iostream>
#include <vector>
#include <vector>
#include <string>


using namespace std;

typedef vector<double> M1;
typedef vector<vector<double> > M2;
typedef vector<vector<vector<double> > > M3;


int main()
{

	int N;		// number of particles
	int navg;	// number of averages
	int seed;	// seed for the random number gen.
	double B;	// magnetic field strength
	double dt;	// time step for the integration
	double tsamp; // time between calculating averages
	double L;	// box size
	int nprint; // print status every nprint

	
	const string inputName = "input.txt";

	ConfigFile config(inputName);
	N = config.read<int>("N");
	navg = config.read<int>("navg");
	seed = config.read<int>("seed");
	B = config.read<double>("B");
	dt = config.read<double>("dt");
	tsamp = config.read<double>("tsamp");
	L = config.read<double>("L");
	nprint = config.read<int>("nprint");


	double m0 = 0.2;
	double m1 = 0.02;
	double m2 = 0.002;

	// the random number generator
	Ranq2 ranNR0(seed);
	Ranq2 ranNR2(seed);
	Ranq2 ranNR1(seed);
	Ranq2 ranNR(seed);

	// vectors describing the state of the system
	M2 r0(N,M1(3));
	M2 dr0(N,M1(3));
	M2 v0(N,M1(3));
	M2 r1(N,M1(3));
	M2 dr1(N,M1(3));
	M2 v1(N,M1(3));
	M2 r2(N,M1(3));
	M2 dr2(N,M1(3));
	M2 v2(N,M1(3));

	M2 r(N,M1(3));
	M2 dr(N,M1(3));

	init_r(r,L,0,ranNR);
	init_r(r0,L,0,ranNR0);
	init_r(r1,L,0,ranNR1);
	init_r(r2,L,0,ranNR2);

	M3 r0t(navg,M2(N,M1(3)));
	M3 r1t(navg,M2(N,M1(3)));
	M3 r2t(navg,M2(N,M1(3)));
	M3 rt(navg,M2(N,M1(3)));


	// use a seperate random number generator
	// for the init of the v's.
	Ranq2 ranNRv0(seed-10);
	Ranq2 ranNRv1(seed-10);
	Ranq2 ranNRv2(seed-10);
	init_v(v0,0.01,ranNRv0);
	init_v(v1,0.01,ranNRv1);
	init_v(v2,0.01,ranNRv2);


	// the deriv objecs takes care of 
	// a single step of dt
	Deriv_mass deriv0(N,m0,B,L,ranNR0);
	Deriv_mass deriv1(N,m1,B,L,ranNR1);
	Deriv_mass deriv2(N,m2,B,L,ranNR2);
	Deriv_od deriv(N,B,L,ranNR);


	double t = 0;
	double t0 = 0;
	double t1 = 0;
	double t2 = 0;
	vector<double> tvec(navg);
	for(int n=0;n<navg;++n) {
		if(n%nprint==0)	cout << n << '\t' << navg << endl;
		
		integrate_mass(r0,dr0,v0,deriv0,t0,t0+tsamp,dt);
		integrate_mass(r1,dr1,v1,deriv1,t1,t1+tsamp,dt);
		integrate_mass(r2,dr2,v2,deriv2,t2,t2+tsamp,dt);

		integrate_od(r,dr,deriv,t,t+tsamp,dt);


		r0t[n] = r0;
		r1t[n] = r1;
		r2t[n] = r2;
		rt[n] = r;	
		tvec[n] = t;
	}

	// save data



	
	return 0;

}








