/*
 *	DESCRIPTION: 	This file contains function definitions for various types of ghost
 *			fluid method. These methods set the state of all ghost cells
 *			(generally by solving a multi material Riemann problem), and also
 *			set the extension_interface_velocity which is used to advance the
 *			level set field.
 *
 *	CITATIONS:	
 *			S Sambasivan, H Udaykumar - "Ghost fluid method for strong shock interactions Part 1: fluid - fluid interfaces" - 2009
 * 			L Xu, C Feng, T Liu - "Practical techniques in ghost fluid method for compressible multi-medium flows" - 2016
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
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
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
				
				p1 = state1.eos->p(ghoststate1.CV(i,all));
				u1 = ghoststate1.get_u(i);
				rho1 = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, p1);
				ghoststate1.CV(i,all) = conserved_variables(rho1,u1,p1,state1.eos);
				
				p2 = state2.eos->p(ghoststate2.CV(i+1,all));
				u2 = ghoststate2.get_u(i+1);
				rho2 = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, p2);
				ghoststate2.CV(i+1,all) = conserved_variables(rho2,u2,p2,state2.eos);
				
				assert(is_state_physical(ghoststate1.CV(i,all), state1.eos));
				assert(is_state_physical(ghoststate1.CV(i+1,all), state1.eos));
				assert(is_state_physical(ghoststate2.CV(i,all), state2.eos));
				assert(is_state_physical(ghoststate2.CV(i+1,all), state2.eos));
				
				
				// Set interface advection velocity as real velocity
				
				extension_interface_velocity(i) = state1.get_u(i);
				extension_interface_velocity(i+1) = state2.get_u(i+1);
			}
			else
			{
				/* 
				 *	Set ghost fluid 1 in cell i, and ghost fluid 2 in cell i+1.
				 *	The fluid 1 entropy is extrapolated from cell i+2 into cell i and i+1
				 * 	The fluid 2 entropy is extrapolated from cell i-1 into cell i and i+1
				 */

				double u2 = state1.get_u(i+1);
				double p2 = state1.eos->p(state1.CV(i+1,all));
				
				double u1 = state2.get_u(i);
				double p1 = state2.eos->p(state2.CV(i,all));
				
				double ref_rho1, ref_p1, ref_rho2, ref_p2;
				
				if (fiii <= 0.0)
				{
					ref_rho1 = state1.CV(i+2,0);
					ref_p1 = state1.eos->p(state1.CV(i+2,all));
				}
				else
				{
					ref_rho1 = state1.CV(i+1,0);
					ref_p1 = state1.eos->p(state1.CV(i+1,all));
				}
				
				if (fm > 0.0)
				{
					ref_rho2 = state2.CV(i-1,0);
					ref_p2 = state2.eos->p(state2.CV(i-1,all));
				}
				else
				{
					ref_rho2 = state2.CV(i,0);
					ref_p2 = state2.eos->p(state2.CV(i,all));
				}
				
				double rho2 = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, p2);
				double rho1 = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, p1);
				
				
				// Set new ghost states

				ghoststate2.CV(i+1,all) = conserved_variables(rho2,u2,p2,state2.eos);
				ghoststate1.CV(i,all) = conserved_variables(rho1,u1,p1,state1.eos);
				
				
				// Set real fluid densities using extrapolated entropy (isobaric fix)
				
				ghoststate2.CV(i,0) = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, state2.eos->p(state2.CV(i,all)));
				ghoststate1.CV(i+1,0) = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, state1.eos->p(state1.CV(i+1,all)));
				
				p1 = state1.eos->p(ghoststate1.CV(i+1,all));
				u1 = ghoststate1.get_u(i+1);
				rho1 = state1.eos->rho_constant_entropy(ref_p1, ref_rho1, p1);
				ghoststate1.CV(i+1,all) = conserved_variables(rho1,u1,p1,state1.eos);
				
				p2 = state2.eos->p(ghoststate2.CV(i,all));
				u2 = ghoststate2.get_u(i);
				rho2 = state2.eos->rho_constant_entropy(ref_p2, ref_rho2, p2);
				ghoststate2.CV(i,all) = conserved_variables(rho2,u2,p2,state2.eos);
				
				assert(is_state_physical(ghoststate1.CV(i,all), state1.eos));
				assert(is_state_physical(ghoststate1.CV(i+1,all), state1.eos));
				assert(is_state_physical(ghoststate2.CV(i,all), state2.eos));
				assert(is_state_physical(ghoststate2.CV(i+1,all), state2.eos));
				
				
				// Set interface advection velocity as real velocity
				
				extension_interface_velocity(i) = state2.get_u(i);
				extension_interface_velocity(i+1) = state1.get_u(i+1);
			}
			
			
		}
	}

		
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}



newR_GFM :: newR_GFM (arrayinfo array)
:
	GFM_base	(array)
{}



void newR_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
)
{
	/*
	 *	A slight modification of the RGFM where the ghost states are set as
	 *	a convex combination of states.
	 */
	
	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;
	double p_star, u_star, rho_star_L, rho_star_R, c_star_L, c_star_R, k_L, k_R, c_L, c_R, xL, xR;

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
				
				c_star_L = state1.eos->a(rho_star_L, p_star);
				c_star_R = state2.eos->a(rho_star_R, p_star);
				c_L = state1.eos->a(state1.CV(i,all));
				c_R = state2.eos->a(state2.CV(i+1,all));
				
				xL = (0.5 * (u_star + (state1.CV(i,1)/state1.CV(i,0))) - 0.5 * (c_star_L + c_L)) * dt;
				xR = std::max((u_star + c_star_R) * dt,  ((state2.CV(i+1,1)/state2.CV(i+1,0)) + c_R ) * dt);
				
				k_L = - std::min(0.0, xL) / state1.array.dx;
				k_R = std::max(0.0, xR) / state1.array.dx;
				
				assert(k_L <= 1.0);
				assert(k_R <= 1.0);
				
				// Prepare vectors of conserved variables for integration
				
				blitz::Array<double,1> U_star_L (conserved_variables(rho_star_L, u_star, p_star, state1.eos));
				blitz::Array<double,1> U_star_R (conserved_variables(rho_star_R, u_star, p_star, state2.eos));
				
				blitz::Array<double,1> new_U1 (3);
				blitz::Array<double,1> new_U2 (3);
				
				new_U1 = k_L * U_star_L + (1.0 - k_L) * ghoststate1.CV(i,all);
				new_U2 = k_R * U_star_R + (1.0 - k_R) * ghoststate2.CV(i+1,all);
				
				std::cout << "Using k_L = " << k_L << std::endl;
				std::cout << "Using k_R = " << k_R << std::endl;
				
				//~ prim_star_L(0) = rho_star_L;
				//~ prim_star_L(1) = u_star;
				//~ prim_star_L(2) = p_star;
				//~ 
				//~ prim_star_R(0) = rho_star_R;
				//~ prim_star_R(1) = u_star;
				//~ prim_star_R(2) = p_star;
				//~ 
				//~ blitz::Array<double,1> W_L (3);
				//~ blitz::Array<double,1> W_R (3);
				//~ 
				//~ W_L(0) = ghoststate1.CV(i,0);
				//~ W_L(1) = ghoststate1.CV(i,1) / ghoststate1.CV(i,0);
				//~ W_L(2) = state1.eos->p(ghoststate1.CV(i,all));
				//~ 
				//~ W_R(0) = ghoststate2.CV(i+1,0);
				//~ W_R(1) = ghoststate2.CV(i+1,1) / ghoststate2.CV(i+1,0);
				//~ W_R(2) = state2.eos->p(ghoststate2.CV(i+1,all));
				
				// Set ghost state as convex combination of states

				ghoststate1.CV(i+1,all) = U_star_L;
				ghoststate1.CV(i,all) = new_U1;
				ghoststate2.CV(i,all) = new_U2;
				ghoststate2.CV(i+1,all) = new_U2;
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
					
				c_star_L = state2.eos->a(rho_star_L, p_star);
				c_star_R = state1.eos->a(rho_star_R, p_star);
				
				k_L = - std::min(0.0, (u_star - c_star_L) * dt) / state1.array.dx;
				k_R = std::max(0.0, (u_star + c_star_R) * dt) / state1.array.dx;
				
				assert(k_L <= 1.0);
				assert(k_R <= 1.0);
				
				// Prepare vectors of conserved variables for integration
				
				blitz::Array<double,1> U_star_L (conserved_variables(rho_star_L, u_star, p_star, state2.eos));
				blitz::Array<double,1> U_star_R (conserved_variables(rho_star_R, u_star, p_star, state1.eos));
				
				blitz::Array<double,1> new_U1 (3);
				blitz::Array<double,1> new_U2 (3);
				
				new_U2 = k_L * U_star_L + (1.0 - k_L) * ghoststate2.CV(i,all);
				new_U1 = k_R * U_star_R + (1.0 - k_R) * ghoststate1.CV(i+1,all);

				ghoststate2.CV(i+1,all) = new_U2;
				ghoststate2.CV(i,all) = new_U2;
				ghoststate1.CV(i,all) = new_U1;
				ghoststate1.CV(i+1,all) = new_U1;
			}
			
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
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
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
)
{
	/*
	 *	Implementation of the Riemann ghost fluid method of Sambasivan and Udaykumar.
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
				ghoststate1.CV(i,all) = conserved_variables(rho_star_L, u_star, p_star, state1.eos);
				ghoststate2.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state2.eos);
				ghoststate2.CV(i+1,all) = conserved_variables(rho_star_R, u_star, p_star, state2.eos);
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
				ghoststate2.CV(i,all) = conserved_variables(rho_star_L, u_star, p_star, state2.eos);
				ghoststate1.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state1.eos);
				ghoststate1.CV(i+1,all) = conserved_variables(rho_star_R, u_star, p_star, state1.eos);
			}
			
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
		}
	}
	
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}




M_GFM :: M_GFM (arrayinfo array)
:
	GFM_base	(array)
{}



void M_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
)
{
	/*
	 *	Implementation of the modified ghost fluid method of 
	 * 	Hu and Khoo. A mixed Riemann problem is solved across
	 * 	the interface, and the ghost state is set to the left 
	 * 	star state, while the entropy is extrapolated into
	 * 	the adjacent real cell.
	 */
	
	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;
	double p_star, u_star, rho_star_L, rho_star_R, newrho, newp, newu;

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
				
				
				// Set density of real fluid by extrapolating entropy from mixed Riemann problem state
				
				newp = state1.eos->p(ghoststate1.CV(i,all));
				newu = ghoststate1.get_u(i);
				newrho = state1.eos->rho_constant_entropy(p_star, rho_star_L, newp);
				ghoststate1.CV(i,all) = conserved_variables(newrho, newu, newp, state1.eos);
				
				newp = state2.eos->p(ghoststate2.CV(i+1,all));
				newu = ghoststate2.get_u(i+1);
				newrho = state2.eos->rho_constant_entropy(p_star, rho_star_R, newp);
				ghoststate2.CV(i+1,all) = conserved_variables(newrho, newu, newp, state2.eos);
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
				
				
				// Set density of real fluid by extrapolating entropy from mixed Riemann problem state
				
				newp = state2.eos->p(ghoststate2.CV(i,all));
				newu = ghoststate2.get_u(i);
				newrho = state2.eos->rho_constant_entropy(p_star, rho_star_L, newp);
				ghoststate2.CV(i,all) = conserved_variables(newrho, newu, newp, state2.eos);
				
				newp = state1.eos->p(ghoststate1.CV(i+1,all));
				newu = ghoststate1.get_u(i+1);
				newrho = state1.eos->rho_constant_entropy(p_star, rho_star_R, newp);
				ghoststate1.CV(i+1,all) = conserved_variables(newrho, newu, newp, state1.eos);
			}
			
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
			
			assert(is_state_physical(ghoststate1.CV(i,all), state1.eos));
			assert(is_state_physical(ghoststate1.CV(i+1,all), state1.eos));
			assert(is_state_physical(ghoststate2.CV(i,all), state2.eos));
			assert(is_state_physical(ghoststate2.CV(i+1,all), state2.eos));
		}
	}
	
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}





P_GFM :: P_GFM (arrayinfo array)
:
	GFM_base	(array)
{}



void P_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
)
{
	/*
	 *	Implementation of the practical ghost fluid method of 
	 * 	Xu and Feng. The interface velocity is found by
	 * 	solving a mixed Riemann problem across the interface, 
	 * 	and then the ghost cells are set using the moving wall
	 * 	boundary conditions. The entropy from the mixed Riemann
	 * 	problem is extrapolated to all ghost cells as well as to
	 * 	the outermost real cell, following the PGFM+A algorithm.
	 */
	
	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	double p_star, u_star, rho_star_L, rho_star_R, newu, newrho, newp;
	
	
	/* 
	 * First set the real fluid state of freshly-cleared cells. This is
	 * done by taking the solution of a mixed Riemann problem between the
	 * interfacial states at time level n-1.
	 */
	
	for (int i=state1.array.numGC; i<state1.array.numGC + state1.array.length-1; i++)
	{
		double phim = ls.linear_interpolation(state1.array.cellcentre_coord(i-1));
		double phi = ls.linear_interpolation(state1.array.cellcentre_coord(i));
		double prevphi = ls_prev.linear_interpolation(state1.array.cellcentre_coord(i));
		double phii = ls.linear_interpolation(state1.array.cellcentre_coord(i+1));
		
		if (phi * prevphi <= 0.0)
		{
			if (phi <= 0.0)
			{
				assert(phim <= 0.0 || phii <= 0.0);
				
				if (phim <= 0.0)
				{
					// Assume that phii is positive and compute mixed Riemann problem between i-1 and i+1
					
					assert(phii > 0.0);
					
					RS->solve_rp_forinterfaceboundary(

						ghoststate1.CV(i-1,all),
						ghoststate2.CV(i+1,all),
						p_star,
						u_star,
						rho_star_L,
						rho_star_R,
						state1.eos,
						state2.eos);
					
					state1.CV(i,all) = conserved_variables(rho_star_L, u_star, p_star, state1.eos);
				}
				else
				{
					assert(phii <= 0.0);
					
					RS->solve_rp_forinterfaceboundary(

						ghoststate2.CV(i-1,all),
						ghoststate1.CV(i+1,all),
						p_star,
						u_star,
						rho_star_L,
						rho_star_R,
						state1.eos,
						state2.eos);
					
					state1.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state1.eos);
				}
				
				
			}
			else
			{
				assert(phim > 0.0 || phii > 0.0);
				
				if (phim > 0.0)
				{
					// Assume that phii is negative and compute mixed Riemann problem between i-1 and i+1
					
					assert(phii <= 0.0);
					
					RS->solve_rp_forinterfaceboundary(

						ghoststate2.CV(i-1,all),
						ghoststate1.CV(i+1,all),
						p_star,
						u_star,
						rho_star_L,
						rho_star_R,
						state1.eos,
						state2.eos);
					
					state2.CV(i,all) = conserved_variables(rho_star_L, u_star, p_star, state2.eos);
				}
				else
				{
					assert(phii > 0.0);
					
					RS->solve_rp_forinterfaceboundary(
	
						ghoststate1.CV(i-1,all),
						ghoststate2.CV(i+1,all),
						p_star,
						u_star,
						rho_star_L,
						rho_star_R,
						state1.eos,
						state2.eos);
					
					state2.CV(i,all) = conserved_variables(rho_star_R, u_star, p_star, state2.eos);
				}
			}
		}
	}
	
	
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;
	

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
				
				
				// Set fluid 1 ghost states in cells i+1 and i+2

				ghoststate1.CV(i+1,all) = conserved_variables(	state1.eos->rho_constant_entropy(p_star, rho_star_L, state1.eos->p(state1.CV(i,all))), 
										2.0*u_star - state1.get_u(i), 
										state1.eos->p(state1.CV(i,all)), 
										state1.eos);
										
				ghoststate1.CV(i+2,all) = conserved_variables(	state1.eos->rho_constant_entropy(p_star, rho_star_L, state1.eos->p(state1.CV(i-1,all))), 
										2.0*u_star - state1.get_u(i-1), 
										state1.eos->p(state1.CV(i-1,all)), 
										state1.eos);
				
				
				// Set fluid 2 ghost states in cells i and i-1
				
				ghoststate2.CV(i,all) = conserved_variables(	state2.eos->rho_constant_entropy(p_star, rho_star_R, state2.eos->p(state2.CV(i+1,all))), 
										2.0*u_star - state2.get_u(i+1), 
										state2.eos->p(state2.CV(i+1,all)), 
										state2.eos);
										
				ghoststate2.CV(i-1,all) = conserved_variables(	state2.eos->rho_constant_entropy(p_star, rho_star_R, state2.eos->p(state2.CV(i+2,all))), 
										2.0*u_star - state2.get_u(i+2), 
										state2.eos->p(state2.CV(i+2,all)), 
										state2.eos);
				
				
				// Set density of real fluid by extrapolating entropy from mixed Riemann problem state
				
				newp = state1.eos->p(ghoststate1.CV(i,all));
				newu = ghoststate1.get_u(i);
				newrho = state1.eos->rho_constant_entropy(p_star, rho_star_L, newp);
				ghoststate1.CV(i,all) = conserved_variables(newrho, newu, newp, state1.eos);
				
				newp = state2.eos->p(ghoststate2.CV(i+1,all));
				newu = ghoststate2.get_u(i+1);
				newrho = state2.eos->rho_constant_entropy(p_star, rho_star_R, newp);
				ghoststate2.CV(i+1,all) = conserved_variables(newrho, newu, newp, state2.eos);
				
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
					

				// Set fluid 2 ghost states in cells i+1 and i+2

				ghoststate2.CV(i+1,all) = conserved_variables(	state2.eos->rho_constant_entropy(p_star, rho_star_L, state2.eos->p(state2.CV(i,all))), 
										2.0*u_star - state2.get_u(i), 
										state2.eos->p(state2.CV(i,all)), 
										state2.eos);
										
				ghoststate2.CV(i+2,all) = conserved_variables(	state2.eos->rho_constant_entropy(p_star, rho_star_L, state2.eos->p(state2.CV(i-1,all))), 
										2.0*u_star - state2.get_u(i-1), 
										state2.eos->p(state2.CV(i-1,all)), 
										state2.eos);
				
				
				// Set fluid 1 ghost states in cells i and i-1
				
				ghoststate1.CV(i,all) = conserved_variables(	state1.eos->rho_constant_entropy(p_star, rho_star_R, state1.eos->p(state1.CV(i+1,all))), 
										2.0*u_star - state1.get_u(i+1), 
										state1.eos->p(state1.CV(i+1,all)), 
										state1.eos);
										
				ghoststate1.CV(i-1,all) = conserved_variables(	state1.eos->rho_constant_entropy(p_star, rho_star_R, state1.eos->p(state1.CV(i+2,all))), 
										2.0*u_star - state1.get_u(i+2), 
										state1.eos->p(state1.CV(i+2,all)), 
										state1.eos);
				
				
				// Set density of real fluid by extrapolating entropy from mixed Riemann problem state
				
				ghoststate1.CV(i+1,0) = state1.eos->rho_constant_entropy(p_star, rho_star_R, state1.eos->p(state1.CV(i+1,all)));
				ghoststate2.CV(i,0) = state2.eos->rho_constant_entropy(p_star, rho_star_L, state2.eos->p(state2.CV(i,all)));
				
				newp = state2.eos->p(ghoststate2.CV(i,all));
				newu = ghoststate2.get_u(i);
				newrho = state2.eos->rho_constant_entropy(p_star, rho_star_L, newp);
				ghoststate2.CV(i,all) = conserved_variables(newrho, newu, newp, state2.eos);
				
				newp = state1.eos->p(ghoststate1.CV(i+1,all));
				newu = ghoststate1.get_u(i+1);
				newrho = state1.eos->rho_constant_entropy(p_star, rho_star_R, newp);
				ghoststate1.CV(i+1,all) = conserved_variables(newrho, newu, newp, state1.eos);
			}
			
			extension_interface_velocity(i-1) = u_star;
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
			extension_interface_velocity(i+2) = u_star;
			
			assert(is_state_physical(ghoststate1.CV(i,all), state1.eos));
			assert(is_state_physical(ghoststate1.CV(i+1,all), state1.eos));
			assert(is_state_physical(ghoststate2.CV(i,all), state2.eos));
			assert(is_state_physical(ghoststate2.CV(i+1,all), state2.eos));
		}
	}
	
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity, 2);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}








newmethod1_GFM :: newmethod1_GFM (arrayinfo array)
:
	GFM_base	(array)
{}



void newmethod1_GFM :: set_ghost_cells (

	fluid_state_array& state1,
	fluid_state_array& state2,
	levelset_array& ls, 
	levelset_array& ls_prev,
	std::shared_ptr<multimat_RS_base> RS,
	double dt
)
{
	/*
	 *	A simple new GFM. Solve mixed Riemann problem to get the interface
     *  velocity for updating level set. Set ghost cells equal to real cell
     *  adjacent to interface.
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

				ghoststate1.CV(i+1,all) = ghoststate1.CV(i,all);
				ghoststate2.CV(i,all) = ghoststate2.CV(i+1,all);
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

				ghoststate2.CV(i+1,all) = ghoststate2.CV(i,all);
				ghoststate1.CV(i,all) = ghoststate1.CV(i+1,all);
			}
			
			extension_interface_velocity(i) = u_star;
			extension_interface_velocity(i+1) = u_star;
		}
	}
	
	extension_advection_eqn_1D(ls, ghoststate1, ghoststate2, extension_interface_velocity);

	state1.CV = ghoststate1.CV;
	state2.CV = ghoststate2.CV;
}



