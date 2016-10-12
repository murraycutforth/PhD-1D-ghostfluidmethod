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




double eos_idealgas :: get_Tau (blitz::Array<double,1> state)
{
	return (gamma - 1.0);
}



double eos_idealgas :: get_Psi (blitz::Array<double,1> state)
{
	return (gamma - 1.0)*specific_ie_cv(state);
}




double eos_idealgas :: get_gamma ()
{
	return gamma;
}



double eos_idealgas :: postshock_density (double P_star, double P_K, double rho_K)
{
	// Use R-H conditions to find star state density shock wave

	return rho_K*((P_star/P_K) + ((gamma-1.0)/(gamma+1.0)))/(1.0 + (P_star/P_K)*((gamma-1.0)/(gamma+1.0)));
}



double eos_idealgas :: postrarefaction_density (double P_star, double P_K, double rho_K)
{
	// Use isentropic law to find star state density inside rarefaction wave

	return rho_K*std::pow(P_star/P_K, 1.0/gamma);
}


eos_type eos_idealgas :: get_eos_type ()
{
	return ideal;
}


































eos_tait :: eos_tait (double gamma, double B)
:
	gamma (gamma),
	B (B)
{}




double eos_tait :: a (blitz::Array<double,1> state)
{
	return (gamma/state(0))*(p(state) + B);
}



double eos_tait :: p (blitz::Array<double,1> state)
{
	return (gamma-1)*state(0)*specific_ie_cv(state) - gamma*B;
}


double eos_tait :: E (blitz::Array<double,1> primitives)
{
	return 0.5*primitives(0)*primitives(1)*primitives(1) + primitives(0)*specific_ie_prim(primitives);
}


double eos_tait :: E (double rho, double u, double p)
{
	return 0.5*rho*u*u + rho*specific_ie_prim(rho,u,p);
}


double eos_tait :: rho_pS (double p, double S)
{
	// Unset
	assert(false);
}


double eos_tait :: S_prho (double p, double rho)
{
	// Unset
	assert(false);
}


double eos_tait :: specific_ie_prim (blitz::Array<double,1> primitives)
{
	return (primitives(2) + gamma*B)/((gamma-1.0)*primitives(0));
}


double eos_tait :: specific_ie_prim (double rho, double u, double p)
{
	return (p + gamma*B)/((gamma-1.0)*rho);
}
	

double eos_tait :: get_Tau (blitz::Array<double,1> state)
{
	return (gamma - 1.0);
}


double eos_tait :: get_Psi (blitz::Array<double,1> state)
{
	return (gamma-1.0)*specific_ie_cv(state);
}

double eos_tait :: get_gamma ()
{
	return gamma;
}

	
double eos_tait :: postshock_density (double P_star, double P_K, double rho_K)
{
	return rho_K*(((2.0*gamma*B)/(P_K*(gamma-1.0))) + (P_star/P_K) + ((gamma-1.0)/(gamma+1.0)))/
	(1.0 + ((2.0*gamma*B)/(P_K*(gamma-1.0)))+ (P_star/P_K)*((gamma-1.0)/(gamma+1.0)));
}


double eos_tait :: postrarefaction_density (double P_star, double P_K, double rho_K)
{
	// Not derived yet
	assert(false);
}



eos_type eos_tait :: get_eos_type ()
{
	return tait;
}








