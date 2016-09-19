#include "eos.hpp"



#include <cmath>




eos_base :: eos_base ()
{}





eos_idealgas :: eos_idealgas (double gamma)
:
	eos_base (),
	gamma (gamma)
{}







double specific_ie_cv (blitz::Array<double,1> state)
{
	return (state(2)/state(0)) - 0.5*(state(1)/state(0))*(state(1)/state(0));
}






bool is_state_physical (blitz::Array<double,1> state)
{
	return (state(0) > 0.0) && (specific_ie_cv(state) > 0.0);
}




double eos_idealgas :: a (blitz::Array<double,1> state)
{
	return sqrt(gamma * p(state) / state(0));
}





double eos_idealgas :: p (blitz::Array<double,1> state)
{
	return (gamma - 1)*state(0)*specific_ie_cv(state);
}





double eos_idealgas :: E (blitz::Array<double,1> primitives)
{
	return 0.5*primitives(0)*primitives(1)*primitives(1) + primitives(0)*specific_ie_prim(primitives);
}




double eos_idealgas :: E (double rho, double u, double p)
{
	return 0.5*rho*u*u + rho*specific_ie_prim(rho,u,p);
}




double eos_idealgas :: rho_pS (double p, double S)
{
	return std::pow(p/S, 1.0/gamma);
}



double eos_idealgas :: S_prho (double p, double rho)
{
	return p/std::pow(rho,gamma);
}






double eos_idealgas :: specific_ie_prim (blitz::Array<double,1> primitives)
{
	return primitives(2)/((gamma - 1.0)*primitives(0));
}




double eos_idealgas :: specific_ie_prim (double rho, double u, double p)
{
	return p/((gamma - 1.0)*rho);
}




double eos_idealgas :: get_gamma ()
{
	return gamma;
}



eos_type eos_idealgas :: get_eos_type ()
{
	return ideal;
}
