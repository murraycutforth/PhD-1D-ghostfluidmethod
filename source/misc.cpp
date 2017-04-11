/*
 *	DESCRIPTION:	This file contains various miscellaneous functions which are useful
 *			but don't fit in anywhere else.
 */


#include "misc.hpp"
#include <cmath>


blitz::Array<double,1> euler_flux (
	
	double rho, 
	double u, 
	double P, 
	double E
)
{
	/*
	 *	Compute Euler flux given value of all density, velocity, pressure and total energy
	 */

	blitz::Array<double,1> flux (3);

	flux(0) = rho*u;
	flux(1) = rho*u*u + P;
	flux(2) = u*(E+P);

	return flux;
}
	

blitz::Array<double,1> euler_flux (
	
	blitz::Array<double,1> cv, 
	std::shared_ptr<eos_base> eos
)
{
	/*
	 *	Compute Euler flux given conserved variables and fluid EOS
	 */

	blitz::Array<double,1> flux (3);

	double P = eos->p(cv);
	double u = cv(1)/cv(0); 

	flux(0) = cv(1);
	flux(1) = cv(1)*u + P;
	flux(2) = u*(cv(2)+P);

	return flux;
}


blitz::Array<double,1> conserved_variables (

	double rho, 
	double u, 
	double p, 
	std::shared_ptr<eos_base> eos
)
{
	/*
	 *	Return vector of conserved variables given primitives
	 */

	blitz::Array<double,1> CV (3);

	CV(0) = rho;
	CV(1) = rho*u;
	CV(2) = eos->E(rho,u,p);

	return CV;
}


double specific_ie_cv (blitz::Array<double,1> state)
{
	return (state(2)/state(0)) - 0.5*(state(1)/state(0))*(state(1)/state(0));
}


bool is_state_physical (blitz::Array<double,1> state, std::shared_ptr<eos_base> eos)
{
	return (state(0) > 0.0) && (eos->p(state) > 0.0);
}


double gaussian_function (double A, double mu, double sigma, double x)
{
	return A*std::exp(-((x-mu)*(x-mu))/(2.0*sigma*sigma));
}


bool cell_local_to_interface (int i, arrayinfo& array, levelset_array& ls)
{
	double phi = ls.linear_interpolation(array.cellcentre_coord(i));
	double phiR = ls.linear_interpolation(array.cellcentre_coord(i+1));
	double phiL = ls.linear_interpolation(array.cellcentre_coord(i-1));

	if (phi*phiR <= 0.0 || phi*phiL <= 0.0) return true;
	else return false;
}
