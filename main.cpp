
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "vecmanip.h"
#include "init_random.h"
#include "deriv_mass.h"
#include "deriv_od.h"
#include "integrate.h"

#include <iostream>
#include <iomanip>
#include <fstream>
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
	double m1 = 0.2;
	double m2 = 0.2;

	// the random number generator
	Ranq2 ranNR0(seed);
	Ranq2 ranNR1(seed);
	Ranq2 ranNR2(seed);
	//Ranq2 ranNR(seed);

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

	//M2 r(N,M1(3));
	//M2 dr(N,M1(3));

	//init_r(r,L,0,ranNR);
	init_r(r0,L,0,ranNR0);
	init_r(r1,L,0,ranNR1);
	init_r(r2,L,0,ranNR2);

	M3 r0t(navg,M2(N,M1(3)));
	M3 r1t(navg,M2(N,M1(3)));
	M3 r2t(navg,M2(N,M1(3)));
	//M3 rt(navg,M2(N,M1(3)));


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
	//Deriv_od deriv(N,B,L,ranNR);


	double t = 0;
	double t0 = 0;
	double t1 = 0;
	double t2 = 0;
	vector<double> tvec(navg);
	for(int n=0;n<navg;++n) {
		//if(n%nprint==0)	cout << n << '\t' << navg << endl;
		integrate_mass(r0,dr0,v0,deriv0,t0,t0+tsamp,dt);
		integrate_mass(r1,dr1,v1,deriv1,t1,t1+tsamp,dt);
		integrate_mass(r2,dr2,v2,deriv2,t2,t2+tsamp,dt);
	
		//integrate_od(r,dr,deriv,t,t+tsamp,dt);
		r0t[n] = r0;
		r1t[n] = r1;
		r2t[n] = r2;
		//rt[n] = r;	
		tvec[n] = t0;
	}

	// save data
	ofstream t_out;
	t_out.open("t.dat");

	ofstream x0_out;
	x0_out.open("x0.dat");
	ofstream y0_out;
	y0_out.open("y0.dat");
	ofstream z0_out;
	z0_out.open("z0.dat");

	ofstream x1_out;
	x1_out.open("x1.dat");
	ofstream y1_out;
	y1_out.open("y1.dat");
	ofstream z1_out;
	z1_out.open("z1.dat");

	ofstream x2_out;
	x2_out.open("x2.dat");
	ofstream y2_out;
	y2_out.open("y2.dat");
	ofstream z2_out;
	z2_out.open("z2.dat");

	//ofstream x_out;
	//x_out.open("x.dat");
	//ofstream y_out;
	//y_out.open("y.dat");
	//ofstream z_out;
	//z_out.open("z.dat");

	char sep1 = '\n';
	char sep2 = ',';
	for( int n=0;n<navg;++n) {
		t_out << tvec[n];
		for( int i=0;i<N;++i) {
			x0_out << r0t[n][i][0];
			y0_out << r0t[n][i][1];
			z0_out << r0t[n][i][2];

			x1_out << r1t[n][i][0];
			y1_out << r1t[n][i][1];
			z1_out << r1t[n][i][2];

	
			x2_out << r2t[n][i][0];
			y2_out << r2t[n][i][1];
			z2_out << r2t[n][i][2];

			//x_out << rt[n][i][0];
			//y_out << rt[n][i][1];
			//z_out << rt[n][i][2];
	
			if(i<N-1) {
				x0_out << sep2;
				y0_out << sep2;
				z0_out << sep2;

				x1_out << sep2;
				y1_out << sep2;
				z1_out << sep2;

				x2_out << sep2;
				y2_out << sep2;
				z2_out << sep2;

				//x_out << sep2;
				//y_out << sep2;
				//z_out << sep2;

			}	
		}
		if(n < navg-1) {

			t_out << sep1;

			x0_out << sep1;
			y0_out << sep1;
			z0_out << sep1;

			x1_out << sep1;
			y1_out << sep1;
			z1_out << sep1;

			x2_out << sep1;
			y2_out << sep1;
			z2_out << sep1;

			//x_out << sep1;
			//y_out << sep1;
			//z_out << sep1;
		}
	}

	t_out.close();
	x0_out.close();
	y0_out.close();
	z0_out.close();
	x1_out.close();
	y1_out.close();
	z1_out.close();
	x2_out.close();
	y2_out.close();
	z2_out.close();
	//x_out.close();
	//y_out.close();
	//z_out.close();
	return 0;

}







