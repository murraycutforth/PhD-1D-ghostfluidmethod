/*
 *	DESCRIPTION: 	This file contains function definitions for various types of ghost
 *			fluid method. These methods set the state of all ghost cells
 *			(generally by solving a multi material Riemann problem), and also
 *			set the extension_interface_velocity which is used to advance the
 *			level set field.
 *
 *	CITATIONS:	
 *			S Sambasivan, H Udaykumar - "Ghost fluid method for strong shock interactions Part 1: fluid - fluid interfaces" - 2009
 *
 */

#include "ghost_fluid_method.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "misc.hpp"
#include "extension_advection_eqn.hpp"
#include <cassert>
#include <cmath>


#define all blitz::Range::all()


GFM_base :: GFM_base (arrayinfo array)
:
	extension_interface_velocity (array.length + 2*array.numGC)
{}




Original_GFM :: Original_GFM (arrayinfo array)
:
	GFM_base	(array)
{}


void Original_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	std::shared_ptr<multimat_RS_base> RS
)
{
	/*
	 *	Implementation of the original ghost fluid method. The ghost states are populated
	 *	by copying pressure and velocity from the real cell in that position, and extrapolating
	 *	the entropy across the interface. The isobaric fix is incorporated - extrapolating the
	 * 	entropy from a real cell away from the interface over to the interfacial real state
	 * 	as well as the ghost states.
	 */

	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;

	for (int i=state1.array.numGC; i<state1.array.numGC + state1.array.length-1; i++)
	{
		double fm = ls.linear_interpolation(state1.array.cellcentre_coord(i-1));
		double fi = ls.linear_interpolation(state1.array.cellcentre_coord(i));
		double fii = ls.linear_interpolation(state1.array.cellcentre_coord(i+1));
		double fiii = ls.linear_interpolation(state1.array.cellcentre_coord(i+2));

		if (fi*fii <= 0.0)
		{
			double ustar;

			if (fi <= 0.0)
			{
				/* 
				 *	Set ghost fluid 2 in cell i, and ghost fluid 1 in cell i+1.
				 * 	The fluid 2 entropy is extrapolated from cell i+2 into cell i+1 and i.
				 * 	The fluid 1 entropy is extrapolated from cell i-1 into cell i and i+1.
				 */

				double u2 = state1.get_u(i);
				double p2 = state1.eos->p(state1.CV(i,all));
				
				double u1 = state2.get_u(i+1);
				double p1 = state2.eos->p(state2.CV(i+1,all));
				
				double ref_rho1, ref_p1, ref_rho2, ref_p2;
				
				if (fiii > 0.0)
				{
					ref_rho2 = state2.CV(i+2,0);
					ref_p2 = state2.eos->p(state2.CV(i+2,all));
				}
				else
				{
					ref_rho2 = state2.CV(i+1,0);
					ref_p2 = state2.eos->p(state2.CV(i+1,all));
				}
				
				if (fm <= 0.0)
				{
					ref_rho1 = state1.CV(i-1,0);
					ref_p1 = state1.eos->p(state1.CV(i-1,all));
				}
				else
				{
					ref_rho1 = state1.CV(i,0);
					ref_p1 = state1.eos->p(state1.CV(i,all));
				}
				
				double rho2 = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, p2);
				double rho1 = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, p1);
				
				
				// Set new ghost states

				ghoststate2.CV(i,all) = conserved_variables(rho2,u2,p2,state2.eos);
				ghoststate1.CV(i+1,all) = conserved_variables(rho1,u1,p1,state1.eos);
				
				
				// Set real fluid densities using extrapolated entropy (isobaric fix)
				
				ghoststate1.CV(i,0) = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, state1.eos->p(state1.CV(i,all)));
				ghoststate2.CV(i+1,0) = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, state2.eos->p(state2.CV(i+1,all)));
				
				
				// Set interface advection velocity as real velocity
				
				extension_interface_velocity(i) = state1.get_u(i);
				extension_interface_velocity(i+1) = state2.get_u(i+1);
			}
			else
			{
				// TODO: as above
				/* 
				 *	Set ghost fluid 1 in cell i, and ghost fluid 2 in cell i+1.
				 *	Extrapolate constant entropy across interface to compute ghost density.
				 */

				double u2 = state1.get_u(i+1);
				double p2 = state1.eos->p(state1.CV(i+1,all));
				double u1 = state2.get_u(i);
				double p1 = state2.eos->p(state2.CV(i,all));
				double rho2 = state2.eos->rho_constant_entropy(p1, state2.CV(i,0), p2);
				double rho1 = state1.eos->rho_constant_entropy(p2, state1.CV(i+1,0), p1);

				ghoststate2.CV(i+1,all) = conserved_variables(rho2,u2,p2,state2.eos);
				ghoststate1.CV(i,all) = conserved_variables(rho1,u1,p1,state1.eos);
				ustar = 0.5*(state1.get_u(i+1) + state2.get_u(i));
			}

			// Set extension velocity field in cell i and i+1
			extension_interface_velocity(i) = ustar;
			extension_interface_velocity(i+1) = ustar;
		}
	}

		
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}





R_GFM :: R_GFM (arrayinfo array)
:
	GFM_base	(array)
{}



void R_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	std::shared_ptr<multimat_RS_base> RS
)
{
	/*
	 *	Implementation of the Riemann ghost fluid method of Sambasican and Udaykumar.
	 *	A mixed Riemann problem is solved across the interface, and the resoluting
	 *	star states are used to populate the ghost cells, as well as the real state
	 *	adjacent to the interface.
	 */
	
	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;
	double p_star, u_star, rho_star_L, rho_star_R;

	for (int i=state1.array.numGC; i<state1.array.numGC + state1.array.length-1; i++)
	{
		double phi = ls.linear_interpolation(state1.array.cellcentre_coord(i));
		double phii = ls.linear_interpolation(state1.array.cellcentre_coord(i+1));

		if (phi*phii <= 0.0)
		{
			if (phi <= 0.0)
			{
				/*
				 *	Fluid 1 is on the left of the interface	
				 */
				
				RS->solve_rp_forinterfaceboundary(

					state1.CV(i,all),
					state2.CV(i+1,all),
					p_star,
					u_star,
					rho_star_L,
					rho_star_R,
					state1.eos,
					state2.eos);

				ghoststate1.CV(i+1,all) = conserved_variables(rho_star_L, u_star, p_star, state1.eos);
				ghoststate2.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state2.eos);
			}
			else
			{
				/*
				 *	Fluid2 is on the left of the interface
				 */

				RS->solve_rp_forinterfaceboundary(

					state2.CV(i,all),
					state1.CV(i+1,all),
					p_star,
					u_star,
					rho_star_L,
					rho_star_R,
					state2.eos,
					state1.eos);

				ghoststate2.CV(i+1,all) = conserved_variables(rho_star_L, u_star, p_star, state2.eos);
				ghoststate1.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state1.eos);
			}
			
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
		}
	}
	
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}
