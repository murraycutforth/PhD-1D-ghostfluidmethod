/*
 *	DESCRIPTION:	Various solvers for the single fluid Euler equations. Each solver defines
 *			the 'single_fluid_update' method which advances the conserved variables
 *			in time by +dt.
 *
 *	CITATIONS: 	E Toro - "Riemann solvers and numerical methods for fluid dynamics: A practical introduction" - 1999
 */

#include "flow_solver.hpp"
#include "eos.hpp"
#include "riemann_solver.hpp"


#define all blitz::Range::all()


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





flow_solver_base :: flow_solver_base (std::shared_ptr<singlefluid_RS_base> rs)
:
	rs (rs)
{}




godunov :: godunov (std::shared_ptr<singlefluid_RS_base> rs)
:
	flow_solver_base (rs)
{}




void godunov :: single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt)
{
	/*
	 * Godunov's first order accurate solver
	 */

	assert(oldstate.array == newstate.array);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 1);
	
	double dtodx = dt/oldstate.array.dx;
	blitz::Array<double,1> flux (3);
	newstate.CV = oldstate.CV;


	// Update real cells correctly

	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC + 1; i++)
	{
		rs->solve_rp(oldstate.CV(i-1,all), oldstate.CV(i,all), flux, oldstate.eos);
		
		newstate.CV(i-1, all) -= dtodx*flux;
		newstate.CV(i, all) += dtodx*flux;
	}	
}




MUSCL :: MUSCL (std::shared_ptr<singlefluid_RS_base> rs)
:
	flow_solver_base (rs)
{}


void MUSCL :: single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt)
{
	/*
	 * The second order accurate MUSCL-Hancock solver
	 */

	assert(oldstate.array == newstate.array);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 2);
	
	
	double dtodx = dt/oldstate.array.dx;
	blitz::Array<double,1> flux (3);
	blitz::Array<double,1> L_del_L (3);
	blitz::Array<double,1> L_del_R (3);
	blitz::Array<double,1> R_del_L (3);
	blitz::Array<double,1> R_del_R (3);
	blitz::Array<double,1> slope_L (3);
	blitz::Array<double,1> slope_R (3);
	blitz::Array<double,1> L_BEV_L (3);
	blitz::Array<double,1> L_BEV_R (3);
	blitz::Array<double,1> R_BEV_L (3);
	blitz::Array<double,1> R_BEV_R (3);
	blitz::Array<double,1> L_BEV_R_evolved (3);
	blitz::Array<double,1> R_BEV_L_evolved (3);
	newstate.CV = oldstate.CV;


	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC + 1; i++)
	{
		// Left cell has index i-1, right cell has index i
		
		L_del_L = oldstate.CV(i-1,all) - oldstate.CV(i-2,all);
		L_del_R = oldstate.CV(i,all) - oldstate.CV(i-1,all);
		R_del_L = oldstate.CV(i,all) - oldstate.CV(i-1,all);
		R_del_R = oldstate.CV(i+1,all) - oldstate.CV(i,all);

		slope_L = MUSCL_slope (	0.5, 
					L_del_L,
					L_del_R);
		slope_R = MUSCL_slope (	0.5,
					R_del_L,
					R_del_R);


		// Limit slopes
		
		double beta = 1.0;
		for (int k=0; k<3; k++)
		{
			slope_L(k) = limited_slope(beta,L_del_L(k),L_del_R(k));
			slope_R(k) = limited_slope(beta,R_del_L(k),R_del_R(k));
		}
		

		// Construct boundary extrapolated values for both cells

		L_BEV_L = oldstate.CV(i-1,all) - 0.5*slope_L;
		L_BEV_R = oldstate.CV(i-1,all) + 0.5*slope_L;
		R_BEV_L = oldstate.CV(i,all) - 0.5*slope_R;
		R_BEV_R = oldstate.CV(i,all) + 0.5*slope_R;

		
		// Evolve states by half a time step

		L_BEV_R_evolved = L_BEV_R + 0.5*(dtodx)*(euler_flux(L_BEV_L, oldstate.eos) - euler_flux(L_BEV_R, oldstate.eos));

		R_BEV_L_evolved = R_BEV_L + 0.5*(dtodx)*(euler_flux(R_BEV_L, oldstate.eos) - euler_flux(R_BEV_R, oldstate.eos));


		// If any states are unphysical, revert to zero slope

		if ((!is_state_physical(L_BEV_R_evolved)) || (!is_state_physical(R_BEV_L_evolved)))
		{
			L_BEV_R_evolved = oldstate.CV(i-1,all);
			R_BEV_L_evolved = oldstate.CV(i,all);
		}


		// Update using conventional Riemann problem solution

		rs->solve_rp(L_BEV_R_evolved, R_BEV_L_evolved, flux, oldstate.eos);
		
		newstate.CV(i-1, all) -= dtodx*flux;
		newstate.CV(i, all) += dtodx*flux;

		
		// If any new states are unphysical, repeat with zero slope

		if ((!is_state_physical(newstate.CV(i-1,all))) || (!is_state_physical(newstate.CV(i,all))))
		{
			newstate.CV(i-1,all) += dtodx*flux;
			newstate.CV(i,all) -= dtodx*flux;
			
			L_BEV_R_evolved = oldstate.CV(i-1,all);
			R_BEV_L_evolved = oldstate.CV(i,all);
		
			rs->solve_rp(L_BEV_R_evolved, R_BEV_L_evolved, flux, oldstate.eos);
			newstate.CV(i-1, all) -= dtodx*flux;
			newstate.CV(i, all) += dtodx*flux;
		}
	}
}



blitz::Array<double,1> MUSCL :: MUSCL_slope (

	double omega, 
	blitz::Array<double,1> delU_L, 
	blitz::Array<double,1> delU_R
)
{
	/*
	 *	Return slope of linear subcell reconstruction of conserved variables
	 */

	blitz::Array<double,1> result (3);	
	result = 0.5*(1.0+omega)*delU_L + 0.5*(1.0*omega)*delU_R;
	return result;
}


double MUSCL :: limited_slope (

	double beta, 
	double del_L, 
	double del_R)
{
	/*
	 *	Limit the magnitude of the MUSCL slope reconstruction to reduce oscillations
	 */

	double result;

	if (del_R > 0.0)
	{
		result = std::max(0.0, std::min(beta*del_L, del_R));
		result = std::max(result, std::min(del_L, beta*del_R));
	}
	else
	{
		result = std::min(0.0, std::max(beta*del_L, del_R));
		result = std::min(result, std::max(del_L, beta*del_R));
	}

	return result;
}

