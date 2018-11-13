
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "vecmanip.h"
#include "init_random.h"
#include "deriv.h"
#include "integrate.h"
#include "density.h"
#include "orientation.h"
#include "flux.h"

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
	int nbin;	// number of bins
	double B;	// magnetic field strength
	double w0;	// number of periods of the magnetic field
	double dt;	// time step for the integration
	double tsamp; // time between calculating averages
	double teq;	// equilibration time
	double L;	// box size
	int nprint; // print status every nprint

	
	const string inputName = "input.txt";

	ConfigFile config(inputName);
	N = config.read<int>("N");
	navg = config.read<int>("navg");
	seed = config.read<int>("seed");
	nbin = config.read<int>("nbin");
	B = config.read<double>("B");
	w0 = config.read<double>("w0");
	dt = config.read<double>("dt");
	tsamp = config.read<double>("tsamp");
	teq = config.read<double>("teq");
	L = config.read<double>("L");
	nprint = config.read<int>("nprint");
	double omega = w0*2*acos(-1.)/L;
	double bs = L/(nbin);


	// density in the xy plane
	M2 rho(nbin,M1(nbin,0.));

	// flux in the xy plane
	// in the x direction
	M2 fx0(nbin,M1(nbin,0.));
	// in the y direction
	M2 fy0(nbin,M1(nbin,0.));

	// flux in the xy plane
	// in the x direction
	M2 fx1(nbin,M1(nbin,0.));
	// in the y direction
	M2 fy1(nbin,M1(nbin,0.));

	// flux in the xy plane
	// in the x direction
	M2 fx2(nbin,M1(nbin,0.));
	// in the y direction
	M2 fy2(nbin,M1(nbin,0.));


	// bins
	M1 bins(nbin,0.);
	for(int i=0;i<nbin;++i)
		bins[i] = (i+.5)*bs;


	// the random number generator
	Ranq2 ranNR(seed);

	// vectors describing the state of the system
	M2 r(N,M1(3));
	M2 dr(N,M1(3));
	init_r(r,L,0.,ranNR);
	
	// the deriv objecs takes care of 
	// a single step of dt
	Deriv deriv(N,B,omega,ranNR);

	// integrates uses the deriv object
	// to increment the system teq in time,
	// in steps of dt

	double t = 0;
	integrate(r,dr,deriv,t,t+teq,dt);


	vector<double> x(navg,0);
	vector<double> y(navg,0);
	vector<double> z(navg,0);
	unsigned int nflux = 0;
	double tf;
	for(int n=0;n<navg;++n) {
		if(n%nprint==0)	cout << n << '\t' << navg << endl;
		x[n] = r[0][0];
		y[n] = r[0][1];
		z[n] = r[0][2];

		tf = t+tsamp;
		while( (t+dt) < tf) { 
			deriv(r,dr,t,dt);

			fluxXY(fx0,fy0,r,dr,dr,bs,L,nbin,0.);	
			fluxXY(fx1,fy1,r,dr,dr,bs,L,nbin,0.5);	
			fluxXY(fx2,fy2,r,dr,dr,bs,L,nbin,1.);	
			nflux += 1;
		}

		// integrate the remaining time
		if(t<tf) 
			deriv(r,dr,t,tf-t); 

		density(r,rho,bs,L,nbin);
		
	}

	// normalize the flux and density
	for(int i=0;i<nbin;++i) {
		for(int j=0;j<nbin;++j){
			fx0[i][j] /= dt*bs*bs*nflux;
			fy0[i][j] /= dt*bs*bs*nflux;
			fx1[i][j] /= dt*bs*bs*nflux;
			fy1[i][j] /= dt*bs*bs*nflux;
			fx2[i][j] /= dt*bs*bs*nflux;
			fy2[i][j] /= dt*bs*bs*nflux;

			rho[i][j] /= navg*N*bs*bs;
		}
	}


	cout << endl << "nflux = " << nflux << endl;
	write_vec("x.dat",x);
	write_vec("y.dat",y);
	write_vec("z.dat",z);
	write_vec("bins.dat",bins);
	write_mat("rho.dat",rho);
	write_mat("fx0.dat",fx0);
	write_mat("fy0.dat",fy0);
	write_mat("fx1.dat",fx1);
	write_mat("fy1.dat",fy1);
	write_mat("fx2.dat",fx2);
	write_mat("fy2.dat",fy2);
	
	return 0;

}








