#ifndef GUARD_integrate_h
#define GUARD_integrate_h

#include <vector>

// increments the state of the system over time T,
// in steps of dt

template<class derivobj>
void integrate_mass(
	std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& v,
	derivobj& deriv,double& ti, double tf, double dt)
{
	
	while( (ti+dt) < tf) { 
		deriv(r,dr,v,ti,dt);
	}

	// integrate the remaining time
	if(ti<tf) 
		deriv(r,dr,v,ti,tf-ti); 
}


template<class derivobj>
void integrate_od(
	std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	derivobj& deriv,double& ti, double tf, double dt)
{
	
	while( (ti+dt) < tf) { 
		deriv(r,dr,ti,dt);
	}

	// integrate the remaining time
	if(ti<tf) 
		deriv(r,dr,ti,tf-ti); 
}



#endif
