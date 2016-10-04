#include "timestep.hpp"
#include "flow_solver.hpp"
#include "ghost_fluid_method.hpp"
#include "levelset_advection.hpp"
#include "eos.hpp"
#include "vfield.hpp"
#include "riemann_solver.hpp"



#include <algorithm>
#include <cmath>



#define all blitz::Range::all()




double compute_dt_serial (double CFL, twofluid_array& states, levelset_array& ls, double T, double t)
{
	// Iterate over real signal speeds using level set
	
	double maxu = 0.0;

	for (int i=states.array.numGC; i<states.array.length + states.array.numGC; i++)
	{
		double xpos = states.array.cellcentre_coord(i);

		if (ls(xpos) <= 0.0)
		{
			maxu = std::max(fabs(states.fluid1(i,1)/states.fluid1(i,0) 
					+ states.eos1->a(states.fluid1(i,all))), 
					maxu);
		}
		else
		{
			maxu = std::max(fabs(states.fluid2(i,1)/states.fluid2(i,0)
					+ states.eos2->a(states.fluid2(i,all))), 
					maxu);
		}
	}

	double dt = CFL * states.array.dx/maxu;

	if (t+dt > T) dt = T - t;
	return dt;
}








double compute_dt_serial_onefluid (double CFL, onefluid_array& state, double T, double t)
{
	double maxu;

	for (int i=0; i<state.array.length; i++)
	{
		maxu = std::max(fabs(state.fluid(i,1)/state.fluid(i,0)) + state.eos->a(state.fluid(i,all)), maxu);
	}

	double dt = CFL*state.array.dx/maxu;

	if (t + dt > T) dt = T - t;

	return dt;
}












void advance_time_level (	double dt, 
				twofluid_array& states, 
				twofluid_array& newstates, 
				levelset_array& ls,
				levelset_array& newls,
				std::shared_ptr<ghost_fluid_method_base> gfm,
				std::shared_ptr<flow_solver_base> fs1,
				std::shared_ptr<flow_solver_base> fs2,
				std::shared_ptr<levelset_advection_base> lsadvection,
				std::shared_ptr<vfield_base> vfield,
				std::shared_ptr<riemann_solver_base> RS)
{
	// Set internal ghost fluid cells

	gfm->set_ghost_cells(states, ls, RS);


	// Apply boundary conditions

	states.apply_BCs();


	// Update single fluid problems separately - no boundary conditions enforced yet

	onefluid_array oldstate1 (states.fluid1, states.array, states.eos1);
	onefluid_array oldstate2 (states.fluid2, states.array, states.eos2);
	onefluid_array newstate1 (newstates.fluid1, newstates.array, newstates.eos1);
	onefluid_array newstate2 (newstates.fluid2, newstates.array, newstates.eos2);

	fs1->single_fluid_update(oldstate1, newstate1, dt);
	fs2->single_fluid_update(oldstate2, newstate2, dt);


	// Update level set function

	lsadvection->advection_operator(ls, newls, vfield, dt);



	// At this point the internal values have been set in newls, newstates
}
