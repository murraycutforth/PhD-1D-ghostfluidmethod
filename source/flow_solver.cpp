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
#include "exact_RS_stiffenedgas.hpp"
#include <random>


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
	
	newstate.apply_BCs();
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
	//~ double TOL = 1e-6;
	//~ double rec_rho, rec_u, rec_p;
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
			UR_L = oldstate.CV(i,all);
			UR_R = oldstate.CV(i,all);
		}
		
		if ((!is_state_physical(UR_L, oldstate.eos)) || (!is_state_physical(UR_R, oldstate.eos)))
		{
			UL_L = oldstate.CV(i-1,all);
			UL_R = oldstate.CV(i-1,all);
			UR_L = oldstate.CV(i,all);
			UR_R = oldstate.CV(i,all);
		}
		
		
		// Now estimate the boundary values at time t+0.5dt
		
		U_L = UL_R + 0.5*dtodx*(euler_flux(UL_L, oldstate.eos) - euler_flux(UL_R, oldstate.eos));
		U_R = UR_L + 0.5*dtodx*(euler_flux(UR_L, oldstate.eos) - euler_flux(UR_R, oldstate.eos));
		
		
		// Use zero slope if evolved pressure or density goes negative
		
		if ((!is_state_physical(U_L, oldstate.eos)) || (!is_state_physical(U_R, oldstate.eos)))
		{
			U_L = oldstate.CV(i-1,all);
			U_R = oldstate.CV(i,all);
		}
		
		
		// Use these as input to a Riemann problem as standard
		 
		rs->solve_rp(U_L, U_R, flux, oldstate.eos);
		newstate.CV(i-1,all) -= dtodx*flux;
		newstate.CV(i,all) += dtodx*flux;
		
		
		// Revert to Godunov method if this flux causes unphysical state
		
		if (i > oldstate.array.numGC) 
		{
			if (! is_state_physical(newstate.CV(i-1,all), oldstate.eos))
			{
				std::cout << "Unphysical state in cell " << i-1 << " where new density = " << newstate.CV(i-1,0) << " new velocity = " << newstate.get_u(i-1) << " and pressure = " << newstate.eos->p(newstate.CV(i-1,all)) << std::endl;
				std::cout << "ULL = " << oldstate.CV(i-2,0) << " " << oldstate.CV(i-2,1) << " " << oldstate.CV(i-2,2) << " " << std::endl;
				std::cout << "UL = " << oldstate.CV(i-1,0) << " " << oldstate.CV(i-1,1) << " " << oldstate.CV(i-1,2) << " " << std::endl;
				std::cout << "UR = " << oldstate.CV(i,0) << " " << oldstate.CV(i,1) << " " << oldstate.CV(i,2) << " " << std::endl;
				assert(is_state_physical(newstate.CV(i-1,all), oldstate.eos));
			}
		}
		if (i==oldstate.array.numGC) FL = flux;
		if (i==oldstate.array.length+oldstate.array.numGC) FR = flux;
	}
	
	newstate.apply_BCs();
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







random_choice_method :: random_choice_method (std::shared_ptr<singlefluid_RS_base> rs)
:
	flow_solver_base (rs)
{}




void random_choice_method :: single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt, blitz::Array<double,1> FL, blitz::Array<double,1> FR)
{
	/*
	 * Standard RCM solver - works only with exact Riemann solver.
	 */

	assert(oldstate.array == newstate.array);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 1);
	assert(oldstate.eos->get_eos_type() == "stiff");
		
	blitz::Array<double,1> flux (3);
	newstate.CV = oldstate.CV;
	
	static std::default_random_engine generator {static_cast<long unsigned int>(time(0))};
	static std::uniform_real_distribution<double> dist (0.0,1.0);
	double sample = dist(generator);
	//~ 
	//~ // Try a van der corput sequence
	//~ static int n = 1;
	//~ const int N = 32;
	//~ const int NO_OF_BITS = 5;
	//~ 
	//~ for (int i=0; i < NO_OF_BITS; i++)
	//~ {
		//~ naive ^= 1 << i;
	//~ }
	//~ 
	//~ double sample = 
	
	std::cout << "Using random number = " << sample << std::endl;
	
	blitz::Array<double,1> Lprimitives (3);
	blitz::Array<double,1> Rprimitives (3);
	blitz::Array<double,1> solnprims (3);
 
	
	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC; i++)
	{
		
		// Randomly pick position in cell to sample solution:
				
		if (sample < 0.5)
		{
			// Solve RP between i-1 and i
			
			Lprimitives(0) = oldstate.CV(i-1, 0);
			Lprimitives(1) = oldstate.CV(i-1, 1)/oldstate.CV(i-1, 0);
			Lprimitives(2) = oldstate.eos->p(oldstate.CV(i-1, all));

			Rprimitives(0) = oldstate.CV(i, 0);  
			Rprimitives(1) = oldstate.CV(i, 1)/oldstate.CV(i, 0);
			Rprimitives(2) = oldstate.eos->p(oldstate.CV(i, all));  
		
			exact_rs_stiffenedgas RS (oldstate.eos->get_gamma(), oldstate.eos->get_gamma(), oldstate.eos->get_Pinf(), oldstate.eos->get_Pinf());
			RS.solve_RP(Lprimitives,Rprimitives);
			solnprims = RS.sample_solution(Lprimitives, Rprimitives, oldstate.array.dx * sample / dt);
		}
		else
		{
			// Solve RP between i and i+1
			
			Lprimitives(0) = oldstate.CV(i, 0);
			Lprimitives(1) = oldstate.CV(i, 1)/oldstate.CV(i, 0);
			Lprimitives(2) = oldstate.eos->p(oldstate.CV(i, all));

			Rprimitives(0) = oldstate.CV(i+1, 0);  
			Rprimitives(1) = oldstate.CV(i+1, 1)/oldstate.CV(i+1, 0);
			Rprimitives(2) = oldstate.eos->p(oldstate.CV(i+1, all));  
		
			exact_rs_stiffenedgas RS (oldstate.eos->get_gamma(), oldstate.eos->get_gamma(), oldstate.eos->get_Pinf(), oldstate.eos->get_Pinf());
			RS.solve_RP(Lprimitives,Rprimitives);
			solnprims = RS.sample_solution(Lprimitives, Rprimitives, (sample - 1.0) * oldstate.array.dx / dt);
		}
		
		
		// Set state in this cell to that specified by solnprims
		
		double E = oldstate.eos->E(solnprims);
		newstate.CV(i,0) = solnprims(0);
		newstate.CV(i,1) = solnprims(0) * solnprims(1);
		newstate.CV(i,2) = E;
	}	
	
	newstate.apply_BCs();
}

