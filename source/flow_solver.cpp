/*
 *	DESCRIPTION:	Various solvers for the single fluid Euler equations. Each solver defines
 *			the 'single_fluid_update' method which advances the conserved variables
 *			in time by +dt.
 *
 *	CITATIONS: 	E Toro - "Riemann solvers and numerical methods for fluid dynamics: A practical introduction" - 1999
 */

#include "flow_solver.hpp"
#include "misc.hpp"
#include "eos.hpp"
#include "riemann_solver.hpp"


#define all blitz::Range::all()


flow_solver_base :: flow_solver_base (std::shared_ptr<singlefluid_RS_base> rs)
:
	rs (rs)
{}




godunov :: godunov (std::shared_ptr<singlefluid_RS_base> rs)
:
	flow_solver_base (rs)
{}




void godunov :: single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt, blitz::Array<double,1> FL, blitz::Array<double,1> FR)
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

		if (i==oldstate.array.numGC) FL = flux;
		if (i==oldstate.array.length+oldstate.array.numGC) FR = flux;
	}	
}




MUSCL :: MUSCL (std::shared_ptr<singlefluid_RS_base> rs)
:
	flow_solver_base (rs)
{}


void MUSCL :: single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt, blitz::Array<double,1> FL, blitz::Array<double,1> FR)
{
	/*
	 * The second order accurate MUSCL-Hancock solver
	 */

	assert(oldstate.array == newstate.array);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 2);
	
	
	double dtodx = dt/oldstate.array.dx;
	double TOL = 1e-6;
	double rec_rho, rec_u, rec_p;
	blitz::Array<double,1> flux (3);
	blitz::Array<double,1> diff_L (3);
	blitz::Array<double,1> diff_C (3);
	blitz::Array<double,1> diff_R (3);
	blitz::Array<double,1> UL_L (3);
	blitz::Array<double,1> UL_R (3);
	blitz::Array<double,1> UR_L (3);
	blitz::Array<double,1> UR_R (3);
	blitz::Array<double,1> U_L (3);
	blitz::Array<double,1> U_R (3);
	newstate.CV = oldstate.CV;


	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC + 1; i++)
	{
		/* 
		 * 	Left cell has index i-1, right cell has index i. This iteration 
		 * 	computes the flux across the boundary between these two cells.
		 */
		
		assert(is_state_physical(oldstate.CV(i-1,all), oldstate.eos));
		assert(is_state_physical(oldstate.CV(i,all), oldstate.eos));
		 
		
		// First compute the four boundary extrapolated states with limited slopes - UL_L, UL_R, UR_L, UR_R
		
		diff_L = oldstate.CV(i-1,all) - oldstate.CV(i-2,all);
		diff_C = oldstate.CV(i,all) - oldstate.CV(i-1,all);
		diff_R = oldstate.CV(i+1,all) - oldstate.CV(i,all);
				
		UL_L = oldstate.CV(i-1,all) - 0.5*limited_slope(diff_L, diff_C);
		UL_R = oldstate.CV(i-1,all) + 0.5*limited_slope(diff_L, diff_C);
		UR_L = oldstate.CV(i,all) - 0.5*limited_slope(diff_C, diff_R);
		UR_R = oldstate.CV(i,all) + 0.5*limited_slope(diff_C, diff_R);
		
		
		// Use zero slope if the reconstructed pressure or density goes negative
		
		if ((!is_state_physical(UL_L, oldstate.eos)) || (!is_state_physical(UL_R, oldstate.eos)))
		{
			UL_L = oldstate.CV(i-1,all);
			UL_R = oldstate.CV(i-1,all);
		}
		
		if ((!is_state_physical(UR_L, oldstate.eos)) || (!is_state_physical(UR_R, oldstate.eos)))
		{
			UR_L = oldstate.CV(i,all);
			UR_R = oldstate.CV(i,all);
		}
		
		
		// Now estimate the boundary values at time t+0.5dt
		
		U_L = UL_R + 0.5*dtodx*(euler_flux(UL_L, oldstate.eos) - euler_flux(UL_R, oldstate.eos));
		U_R = UR_L + 0.5*dtodx*(euler_flux(UR_L, oldstate.eos) - euler_flux(UR_R, oldstate.eos));
		
		
		// Enforce a floor on the evolved pressure and density
		
		rec_rho = U_L(0);
		rec_u = U_L(1)/U_L(0);
		rec_p = oldstate.eos->p(U_L);
		rec_rho = std::max(TOL, rec_rho);
		rec_p = std::max(TOL, rec_p);
		U_L = conserved_variables(rec_rho, rec_u, rec_p, oldstate.eos);
		
		rec_rho = U_R(0);
		rec_u = U_R(1)/U_R(0);
		rec_p = oldstate.eos->p(U_R);
		rec_rho = std::max(TOL, rec_rho);
		rec_p = std::max(TOL, rec_p);
		U_R = conserved_variables(rec_rho, rec_u, rec_p, oldstate.eos);
		
		
		// Use these as input to a Riemann problem as standard
		 
		rs->solve_rp(U_L, U_R, flux, oldstate.eos);
		newstate.CV(i-1,all) -= dtodx*flux;
		newstate.CV(i,all) += dtodx*flux;
		 

		if (i==oldstate.array.numGC) FL = flux;
		if (i==oldstate.array.length+oldstate.array.numGC) FR = flux;
	}
}




blitz::Array<double,1> MUSCL :: limited_slope (

	blitz::Array<double,1> del_L, 
	blitz::Array<double,1> del_R
)
{
	/*
	 *	Limit the magnitude of the MUSCL slope reconstruction to reduce oscillations using minbee
	 */

	blitz::Array<double,1> result (3);
	
	for (int k=0; k<3; k++)
	{
		if (del_R(k) > 0.0)
		{
			result(k) = std::max(0.0, std::min(del_L(k), del_R(k)));
		}	
		else
		{
			result(k) = std::min(0.0, std::max(del_L(k), del_R(k)));
		}	
	}

	return result;
}

