#ifndef GUARD_deriv_mass_magnetic_h
#define GUARD_deriv_mass_magnetic_h

#include <iostream>

#include <cmath>
#include <vector>
#include "vecmanip.h"
#include "box_muller.h"

struct Deriv_mass {
	public:

	Deriv_mass(int NN, double mm, double BB, double LL, Ranq2 ranNRR):
		N(NN),m(mm),B(BB), L(LL),ranNR(ranNRR)
		  {sqrt2 = std::sqrt(2.);
			omega = 2*acos(-1.)/L;}


	// evolves the state of the system to t+dt
	void  operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		double& t, double dt);

	private:
		int N;		// number of particles
		double m;	// mass of the particles
		double B;	
		double L;
		double omega;
public:
		Ranq2 ranNR;	// the ranom number generator
private:
		double sqrt2;

		double Bfield(const std::vector<double>& ri) {
				return B*sin(ri[1]*omega);}

		double BfieldP(const std::vector<double>& ri) {
				return B*omega*cos(ri[1]*omega);}
};



void Deriv_mass::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		double& t, double dt)
{
	double sqrt_dt = std::sqrt(dt);

	double Bri; 

	for(int i=0;i<N;++i) {

		Bri = Bfield(r[i]);

		dr[i][0] = v[i][0]*dt;
		dr[i][1] = v[i][1]*dt;
		dr[i][2] = v[i][2]*dt;
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];


		v[i][0] += (-Bri*v[i][1]*dt - v[i][0]*dt + ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][1] += (Bri*v[i][0]*dt - v[i][1]*dt + ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][2] += (-v[i][2]*dt+ ndist(ranNR)*sqrt_dt*sqrt2)/m;	


	}
	
	t += dt;
}





#endif


