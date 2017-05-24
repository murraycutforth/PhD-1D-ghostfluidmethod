/*
 *	DESCRIPTION:	This file contains definitions of all functions which involve the
 *			equation of state of the material.
 *	
 *	TODO:		(High priority) Implement Tait eos
 */


#include "eos.hpp"
#include "misc.hpp"
#include <cmath>



eos_base :: eos_base ()
{}





eos_idealgas :: eos_idealgas (double gamma)
:
	eos_base (),
	gamma (gamma)
{}



std::string eos_idealgas :: get_eos_type ()
{
	return "ideal";
}


double eos_idealgas :: get_gamma ()
{
	return gamma;
}


double eos_idealgas :: get_Pinf ()
{
	return 0.0;
}


double eos_idealgas :: a (blitz::Array<double,1> state)
{
	/*
	 *	Sound speed from conserved variables
	 */

	return sqrt(gamma * p(state) / state(0));
}


double eos_idealgas :: a (double rho, double p)
{
	/*
	 *	Sound speed from conserved variables
	 */

	return sqrt(gamma * p / rho);
}


double eos_idealgas :: p (blitz::Array<double,1> state)
{
	/*
	 *	Pressure from conserved variables
	 */

	return (gamma - 1)*state(0)*specific_ie_cv(state);
}


double eos_idealgas :: E (blitz::Array<double,1> primitives)
{
	/*
	 *	Total energy from primitive variables
	 */

	return 0.5*primitives(0)*primitives(1)*primitives(1) + primitives(0)*specific_ie_prim(primitives);
}


double eos_idealgas :: E (double rho, double u, double p)
{
	/*
	 *	Total energy from primitive variables
	 */

	return 0.5*rho*u*u + rho*specific_ie_prim(rho,u,p);
}


double eos_idealgas :: specific_ie_prim (blitz::Array<double,1> primitives)
{
	/*
	 *	Specific internal energy as a function of primitives
	 */

	return primitives(2)/((gamma - 1.0)*primitives(0));
}


double eos_idealgas :: specific_ie_prim (double rho, double u, double p)
{
	/*
	 *	Specific internal energy as a function of primitives
	 */

	return p/((gamma - 1.0)*rho);
}



double eos_idealgas :: rho_constant_entropy (double p_old, double rho_old, double p_new)
{
	/*
	 *	Density at pressure = p_new along isentrope from state (p_old, rho_old)
	 */

	double entropy = S(p_old, rho_old);
	return rho(p_new, entropy);
}

double eos_idealgas :: rho (double p, double S)
{
	/*
	 *	Density as a function of pressure and entropy
	 */

	return std::pow(p/S, 1.0/gamma);
}


double eos_idealgas :: S (double p, double rho)
{
	/*
	 *	Entropy as a function of pressure and density
	 */
	
	return p/std::pow(rho,gamma);
}



double eos_idealgas :: get_Tau (blitz::Array<double,1> state)
{
	return (gamma - 1.0);
}


double eos_idealgas :: get_Psi (blitz::Array<double,1> state)
{
	return (gamma - 1.0)*specific_ie_cv(state);
}


double eos_idealgas :: postshock_density (double P_star, double P_K, double rho_K)
{
	/* 
	 *	Use R-H conditions to find star state density behind shock wave
	 */

	return rho_K*((P_star/P_K) + ((gamma-1.0)/(gamma+1.0)))/(1.0 + (P_star/P_K)*((gamma-1.0)/(gamma+1.0)));
}


double eos_idealgas :: postrarefaction_density (double P_star, double P_K, double rho_K)
{
	/* 
	 *	Use isentropic law to find star state density behind rarefaction wave
	 */

	return rho_K*std::pow(P_star/P_K, 1.0/gamma);
}















eos_stiffenedgas :: eos_stiffenedgas (double gamma, double Pinf)
:
	eos_base (),
	gamma (gamma),
	Pinf (Pinf)
{}



std::string eos_stiffenedgas :: get_eos_type ()
{
	return "stiff";
}


double eos_stiffenedgas :: get_gamma ()
{
	return gamma;
}

double eos_stiffenedgas :: get_Pinf ()
{
	return Pinf;
}


double eos_stiffenedgas :: a (blitz::Array<double,1> state)
{
	/*
	 *	Sound speed from conserved variables
	 */

	return sqrt(gamma * (p(state) + Pinf) / state(0));
}


double eos_stiffenedgas :: a (double rho, double p)
{
	/*
	 *	Sound speed from conserved variables
	 */

	return sqrt(gamma * (p + Pinf) / rho);
}


double eos_stiffenedgas :: p (blitz::Array<double,1> state)
{
	/*
	 *	Pressure from conserved variables
	 */

	return (gamma - 1)*state(0)*specific_ie_cv(state) - gamma*Pinf;
}


double eos_stiffenedgas :: E (blitz::Array<double,1> primitives)
{
	/*
	 *	Total energy from primitive variables
	 */

	return 0.5*primitives(0)*primitives(1)*primitives(1) + primitives(0)*specific_ie_prim(primitives);
}


double eos_stiffenedgas :: E (double rho, double u, double p)
{
	/*
	 *	Total energy from primitive variables
	 */

	return 0.5*rho*u*u + rho*specific_ie_prim(rho,u,p);
}


double eos_stiffenedgas :: specific_ie_prim (blitz::Array<double,1> primitives)
{
	/*
	 *	Specific internal energy as a function of primitives
	 */

	return (primitives(2) + gamma*Pinf)/((gamma - 1.0)*primitives(0));
}


double eos_stiffenedgas :: specific_ie_prim (double rho, double u, double p)
{
	/*
	 *	Specific internal energy as a function of primitives
	 */

	return (p + gamma*Pinf)/((gamma - 1.0)*rho);
}



double eos_stiffenedgas :: rho_constant_entropy (double p_old, double rho_old, double p_new)
{
	/*
	 *	Density at pressure = p_new along isentrope from state (p_old, rho_old)
	 */

	double entropy = S(p_old, rho_old);
	return rho(p_new, entropy);
}

double eos_stiffenedgas :: rho (double p, double S)
{
	/*
	 *	Density as a function of pressure and entropy
	 */

	return std::pow((p + Pinf)/S, 1.0/gamma);
}


double eos_stiffenedgas :: S (double p, double rho)
{
	/*
	 *	Entropy as a function of pressure and density
	 */
	
	return (p + Pinf)/std::pow(rho,gamma);
}



// Following equations still need to be extended from the ideal case - only used by the M-HLLC solver

double eos_stiffenedgas :: get_Tau (blitz::Array<double,1> state)
{
	return (gamma - 1.0);
}


double eos_stiffenedgas :: get_Psi (blitz::Array<double,1> state)
{
	return (gamma - 1.0)*specific_ie_cv(state);
}


double eos_stiffenedgas :: postshock_density (double P_star, double P_K, double rho_K)
{
	/* 
	 *	Use R-H conditions to find star state density behind shock wave
	 */

	return rho_K*((P_star/P_K) + ((gamma-1.0)/(gamma+1.0)))/(1.0 + (P_star/P_K)*((gamma-1.0)/(gamma+1.0)));
}


double eos_stiffenedgas :: postrarefaction_density (double P_star, double P_K, double rho_K)
{
	/* 
	 *	Use isentropic law to find star state density behind rarefaction wave
	 */

	return rho_K*std::pow(P_star/P_K, 1.0/gamma);
}





