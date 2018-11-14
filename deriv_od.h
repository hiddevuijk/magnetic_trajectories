#ifndef GUARD_deriv_od_magnetic_h
#define GUARD_deriv_od_magnetic_h

#include <iostream>

#include <cmath>
#include <vector>
#include "vecmanip.h"
#include "box_muller.h"

struct Deriv_od {
	public:

	Deriv_od(int NN, double BB, double LL,  
		 Ranq2 ranNRR):
		N(NN),B(BB), L(LL),ranNR(ranNRR)
		  {sqrt2 = std::sqrt(2.);
			omega = 2*acos(-1.)/L;}


	// evolves the state of the system to t+dt
	void  operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double& t, double dt);

	private:
		int N;		// number of particles
		double B;
		double L;
		double omega;
		Ranq2 ranNR;	// the ranom number generator

		double sqrt2;

		double Bfield(const std::vector<double>& ri) {
					return B*sin(ri[1]*omega);}
		double BfieldP(const std::vector<double>& ri) {
					return B*omega*cos(ri[1]*omega);}

};



void Deriv_od::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double& t, double dt)
{
	double sqrt_dt = std::sqrt(dt);
	double etaX, etaY;
	double Bi,Bpi;
	double D;
	for(int i=0;i<N;++i) {

		Bi = Bfield(r[i]);
		Bpi = BfieldP(r[i]);

		D = 1./(1+Bi*Bi);


		etaX = ndist(ranNR)*sqrt_dt*sqrt2;
		etaY = ndist(ranNR)*sqrt_dt*sqrt2;

		dr[i][0] = D*D*(1-Bi*Bi)*Bpi*dt + D*etaX + Bi*D*etaY; 
		dr[i][1] = -2*D*D*Bi*Bpi*dt + D*etaY - Bi*D*etaX;
		dr[i][2] = ndist(ranNR)*sqrt_dt*sqrt2;	

		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

	}
	
	t += dt;
}





#endif


