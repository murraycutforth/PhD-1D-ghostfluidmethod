/*
 *	DESCRIPTION: 	This file contains function definitions for various types of ghost
 *			fluid method. These methods set the state of all ghost cells
 *			(generally by solving a multi material Riemann problem), and also
 *			set the extension_interface_velocity which is used to advance the
 *			level set field.
 *
 *	CITATIONS:
 *
 */

#include "ghost_fluid_method.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "misc.hpp"
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
	 *	the entropy across the interface.
	 */

	assert(ls.array.numGC >= 1);
	assert(state1.array.numGC >= 1);
	static fluid_state_array ghoststate1 (state1.copy());
	static fluid_state_array ghoststate2 (state2.copy());
	ghoststate1.CV = state1.CV;
	ghoststate2.CV = state2.CV;

	for (int i=state1.array.numGC; i<state1.array.numGC + state1.array.length-1; i++)
	{
		double fi = ls.linear_interpolation(state1.array.cellcentre_coord(i));
		double fii = ls.linear_interpolation(state1.array.cellcentre_coord(i+1));

		if (fi*fii <= 0.0)
		{
			double ustar;

			if (fi <= 0.0)
			{
				/* 
				 *	Set ghost fluid 2 in cell i, and ghost fluid 1 in cell i+1.
				 *	Extrapolate constant entropy across interface to compute ghost density.
				 */

				double u2 = state1.get_u(i);
				double p2 = state1.eos->p(state1.CV(i,all));
				double u1 = state2.get_u(i+1);
				double p1 = state2.eos->p(state2.CV(i+1,all));
				double rho2 = state2.eos->rho_constant_entropy(p1, state2.CV(i+1,0), p2);
				double rho1 = state1.eos->rho_constant_entropy(p2, state1.CV(i,0), p1);

				ghoststate2.CV(i,all) = conserved_variables(rho2,u2,p2,state2.eos);
				ghoststate1.CV(i+1,all) = conserved_variables(rho1,u1,p1,state1.eos);
				ustar = 0.5*(state1.get_u(i) + state2.get_u(i+1));
			}
			else
			{
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

		
	// Advect ghost fluid and extension v states outwards (different grids)





	/* OLD CODE
	
	for (int i=ls.array.numGC; i<ls.array.length+ls.array.numGC; i++)
	{

		// Look for changes in sign of level set field (these are material interfaces)

		if (std::copysign(1.0,ls.phi(i)) != std::copysign(1.0,ls.phi(i+1)))
		{
			double gradphi = (ls.phi(i+1) - ls.phi(i))/ls.array.dx;

			
			// Use linear interpolation to find the zero level set

			double interfacepos = -ls.phi(i)/gradphi + ls.array.cellcentre_coord(i);
			assert(interfacepos >= ls.array.cellcentre_coord(i));
			assert(interfacepos <= ls.array.cellcentre_coord(i+1));


			// Identify the indices of fluid cells to the L and R of the interface
			
			int cellindex = states.array.cellindex(interfacepos);

			int L_index, R_index;

			if (interfacepos > states.array.cellcentre_coord(cellindex))
			{
				L_index = cellindex;
				R_index = cellindex + 1;
			}
			else
			{
				L_index = cellindex - 1;
				R_index = cellindex;
			}

			
			blitz::Array<double,1> fluid1realstate (3);
			blitz::Array<double,1> fluid2realstate (3);

			if (gradphi < 0.0)
			{
				fluid1realstate = states.fluid1(R_index,all);
				fluid2realstate = states.fluid2(L_index,all);
			}
			else
			{
				fluid1realstate = states.fluid1(L_index,all);
				fluid2realstate = states.fluid2(R_index,all);
			}


			// Iterate to the left until the gradient of phi changes sign 
			
			double phi_L, phi_R, currentgradphi;

			for (int k=L_index; k>=states.array.numGC; k--)
			{
				phi_L = ls(states.array.cellcentre_coord(k));
				phi_R = ls(states.array.cellcentre_coord(k+1));
				currentgradphi = (phi_R - phi_L)/states.array.dx;

				if (std::copysign(1.0,currentgradphi) == std::copysign(1.0,gradphi))
				{
					// Set ghost cells at index k

					if (gradphi <= 0.0)
					{
						// Fluid 1 is ghost here

						assert(ls(states.array.cellcentre_coord(k)) >= 0.0);

						
						// Copy pressure and velocity from fluid 2 at this cell

						double p = states.eos2->p(states.fluid2(k,all));
						double u = states.fluid2(k,1)/states.fluid2(k,0);


						// Extrapolate entropy from nearest real fluid 1 cell

						double S = states.eos1->S_prho(	states.eos1->p(fluid1realstate),
										fluid1realstate(0));
						double rho = states.eos1->rho_pS(p,S);


						// Set conserved variables

						states.fluid1(k,0) = rho;
						states.fluid1(k,1) = rho*u;
						states.fluid1(k,2) = states.eos1->E(rho,u,p);
					}
					else
					{
						// Fluid 2 is ghost

						assert(ls(states.array.cellcentre_coord(k)) <= 0.0);

						
						// Copy pressure and velocity from fluid 1 at this cell

						double p = states.eos1->p(states.fluid1(k,all));
						double u = states.fluid1(k,1)/states.fluid1(k,0);


						// Extrapolate entropy from nearest fluid 2 cell

						double S = states.eos2->S_prho(	states.eos2->p(fluid2realstate),
										fluid2realstate(0));
						double rho = states.eos2->rho_pS(p,S);


						// Set conserved variables

						states.fluid2(k,0) = rho;
						states.fluid2(k,1) = rho*u;
						states.fluid2(k,2) = states.eos2->E(rho,u,p);
					}
				}
				else
				{
					break;
				}
			}


			// Now repeat, iterating to the right

			for (int k=R_index; k<states.array.length + states.array.numGC; k++)
			{
				phi_L = ls(states.array.cellcentre_coord(k-1));
				phi_R = ls(states.array.cellcentre_coord(k));
				currentgradphi = (phi_R - phi_L)/states.array.dx;

				if (std::copysign(1.0,currentgradphi) == std::copysign(1.0,gradphi))
				{
					// Set ghost cells at index k

					if (gradphi <= 0.0)
					{
						// Fluid 2 is ghost

						assert(ls(states.array.cellcentre_coord(k)) <= 0.0);

						
						// Copy pressure and velocity from fluid 1 at this cell

						double p = states.eos1->p(states.fluid1(k,all));
						double u = states.fluid1(k,1)/states.fluid1(k,0);


						// Extrapolate entropy from nearest fluid 2 cell

						double S = states.eos2->S_prho(	states.eos2->p(fluid2realstate),
										fluid2realstate(0));
						double rho = states.eos2->rho_pS(p,S);


						// Set conserved variables

						states.fluid2(k,0) = rho;
						states.fluid2(k,1) = rho*u;
						states.fluid2(k,2) = states.eos2->E(rho,u,p);
					}
					else
					{
						// Fluid 1 is ghost here

						assert(ls(states.array.cellcentre_coord(k)) >= 0.0);

						
						// Copy pressure and velocity from fluid 2 at this cell

						double p = states.eos2->p(states.fluid2(k,all));
						double u = states.fluid2(k,1)/states.fluid2(k,0);


						// Extrapolate entropy from nearest real fluid 1 cell

						double S = states.eos1->S_prho(	states.eos1->p(fluid1realstate),
										fluid1realstate(0));
						double rho = states.eos1->rho_pS(p,S);


						// Set conserved variables

						states.fluid1(k,0) = rho;
						states.fluid1(k,1) = rho*u;
						states.fluid1(k,2) = states.eos1->E(rho,u,p);
					}
				}
				else
				{
					break;
				}
			}
		}
	}*/
}

















