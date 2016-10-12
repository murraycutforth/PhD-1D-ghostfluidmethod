#include "ghost_fluid_method.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"



#include <cassert>
#include <cmath>



#define all blitz::Range::all()





ghost_fluid_method_base :: ghost_fluid_method_base (arrayinfo array)
:
	extension_interface_velocity (array.length + 2*array.numGC)
{}








original_GFM :: original_GFM (arrayinfo array)
:
	ghost_fluid_method_base(array)
{}


void original_GFM :: set_ghost_cells (	twofluid_array& states,
					levelset_array& ls,
					std::shared_ptr<riemann_solver_base> rs)
{


	assert(ls.array.numGC >= 1);
	
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
	}
}


















isobaric_fix_GFM :: isobaric_fix_GFM (arrayinfo array)
:
	ghost_fluid_method_base (array)
{}


void isobaric_fix_GFM :: set_ghost_cells (	twofluid_array& states,
					levelset_array& ls,
					std::shared_ptr<riemann_solver_base> rs)
{


	assert(ls.array.numGC >= 1);
	
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


			// Isobaric fix - use points one further away from interface

			L_index--;
			R_index++;

			
			blitz::Array<double,1> fluid1realstate (3);
			blitz::Array<double,1> fluid2realstate (3);

			if (gradphi < 0.0)
			{
				fluid1realstate = states.fluid1(R_index,all);
				fluid2realstate = states.fluid2(L_index,all);


				// Extrapolate entropy to real fluids adjacent to interface

				double S1 = states.eos1->S_prho(	states.eos1->p(fluid1realstate),
										fluid1realstate(0));
				double S2 = states.eos2->S_prho(	states.eos2->p(fluid2realstate),
										fluid2realstate(0));

				states.fluid1(R_index-1,0) = states.eos1->rho_pS(states.eos1->p(states.fluid1(R_index-1,all)),S1);
				states.fluid2(L_index+1,0) = states.eos2->rho_pS(states.eos2->p(states.fluid2(L_index+1,all)),S2);

			}
			else
			{
				fluid1realstate = states.fluid1(L_index,all);
				fluid2realstate = states.fluid2(R_index,all);


				// Extrapolate entropy to real fluids adjacent to interface
				
				double S1 = states.eos1->S_prho(	states.eos1->p(fluid1realstate),
										fluid1realstate(0));
				double S2 = states.eos2->S_prho(	states.eos2->p(fluid2realstate),
										fluid2realstate(0));

				states.fluid1(L_index+1,0) = states.eos1->rho_pS(states.eos1->p(states.fluid1(L_index+1,all)),S1);
				states.fluid2(R_index-1,0) = states.eos2->rho_pS(states.eos2->p(states.fluid2(R_index-1,all)),S2);
			}


			// Iterate to the left until the gradient of phi changes sign 
			
			double phi_L, phi_R, currentgradphi;

			for (int k=L_index + 1; k>=states.array.numGC; k--)
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

			for (int k=R_index - 1; k<states.array.length + states.array.numGC; k++)
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
	}
}















rGFM :: rGFM (arrayinfo array)
:
	ghost_fluid_method_base (array)
{}



void rGFM :: set_ghost_cells (	twofluid_array& states,
				levelset_array& ls,
				std::shared_ptr<riemann_solver_base> rs)
{


	assert(ls.array.numGC >= 1);
	assert(ls.array.length+2*ls.array.numGC == extension_interface_velocity.extent(blitz::firstDim));
	
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
			double rho_IL, rho_IR, p_I, u_I;
			blitz::Array<double,1> fluid1interfacestate (3);
			blitz::Array<double,1> fluid2interfacestate (3);

			if (gradphi <= 0.0)
			{
				fluid1realstate = states.fluid1(R_index,all);
				fluid2realstate = states.fluid2(L_index,all);

				// Solve RP with fluid1 on R

				rs->solve_rp_forinterfaceboundary (	fluid2realstate,
									fluid1realstate,
									p_I,
									u_I,
									rho_IL,
									rho_IR,
									states.eos2,
									states.eos1 );

				fluid1interfacestate(0) = rho_IR;
				fluid1interfacestate(1) = rho_IR*u_I;
				fluid1interfacestate(2) = states.eos1->E(rho_IR, u_I, p_I);
				assert(is_state_physical(fluid1interfacestate));

				fluid2interfacestate(0) = rho_IL;
				fluid2interfacestate(1) = rho_IL*u_I;
				fluid2interfacestate(2) = states.eos2->E(rho_IL, u_I, p_I);
				assert(is_state_physical(fluid2interfacestate));


				// Redefine the real fluid next to the interface

				states.fluid2(L_index,all) = fluid2interfacestate(all);
				states.fluid1(R_index,all) = fluid1interfacestate(all);
			}
			else
			{
				fluid1realstate = states.fluid1(L_index,all);
				fluid2realstate = states.fluid2(R_index,all);

				// Solve RP with fluid1 on L

				rs->solve_rp_forinterfaceboundary (	fluid1realstate,
									fluid2realstate,
									p_I,
									u_I,
									rho_IL,
									rho_IR,
									states.eos1,
									states.eos2 );
				
				fluid1interfacestate(0) = rho_IL;
				fluid1interfacestate(1) = rho_IL*u_I;
				fluid1interfacestate(2) = states.eos1->E(rho_IL, u_I, p_I);
				assert(is_state_physical(fluid1interfacestate));

				fluid2interfacestate(0) = rho_IR;
				fluid2interfacestate(1) = rho_IR*u_I;
				fluid2interfacestate(2) = states.eos2->E(rho_IR, u_I, p_I);
				assert(is_state_physical(fluid2interfacestate));


				// Redefine the real fluid next to the interface

				states.fluid1(L_index,all) = fluid1interfacestate(all);
				states.fluid2(R_index,all) = fluid2interfacestate(all);
			}

			

			// Iterate through extension velocity array - filling with closest velocity

			extension_interface_velocity = -1e100;

			for (int k=i; k<=ls.array.length+2*ls.array.numGC-2; k++)
			{
				double thisgradphi = (ls.phi(k+1)-ls.phi(k))/ls.array.dx;

				if (std::copysign(1.0,gradphi) == std::copysign(1.0,thisgradphi) || fabs(thisgradphi) < 1e-10)
				{
					extension_interface_velocity(k) = u_I;
				}
			}
			for (int k=i-1; k>=0; k--)
			{
				double thisgradphi = (ls.phi(k+1) - ls.phi(k))/ls.array.dx;

				if (std::copysign(1.0,gradphi) == std::copysign(1.0,thisgradphi) || fabs(thisgradphi) < 1e-10)
				{
					extension_interface_velocity(k) = u_I;
				}
			}

			// Finally copy value across to rightmost cell (which is undefined so far)

			extension_interface_velocity(ls.array.length+2*ls.array.numGC-1) = extension_interface_velocity(ls.array.length+2*ls.array.numGC-2);


			// Iterate to the left until the gradient of phi changes sign (setting ghost states)
			
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

						states.fluid1(k,all) = fluid1interfacestate(all);

					}
					else
					{
						// Fluid 2 is ghost

						assert(ls(states.array.cellcentre_coord(k)) <= 0.0);
						
						states.fluid2(k,all) = fluid2interfacestate(all);
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
						
						states.fluid2(k,all) = fluid2interfacestate(all);
					}
					else
					{
						// Fluid 1 is ghost here

						assert(ls(states.array.cellcentre_coord(k)) >= 0.0);

						states.fluid1(k,all) = fluid1interfacestate(all);
					}
				}
				else
				{
					break;
				}
			}
		}
	}


	for (int k=0; k<=ls.array.length+2*ls.array.numGC-1; k++)
	{
		assert(extension_interface_velocity(k) != -1e100);
	}
}




