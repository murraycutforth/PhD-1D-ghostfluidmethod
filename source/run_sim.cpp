/*
 *	DESCRIPTION:	File defines classes which implement various simulation modes. Each definition
 *			of the run_sim member function performs all steps required to solve an
 *			initial-boundary value problem for the Euler equations, and output results
 *			to file.
 *
 *	TODO:		(High priority) Refactor the twofluid_sim code
 *
 */


#include "run_sim.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "misc.hpp"
#include "construct_initialise.hpp"
#include "error.hpp"
#include <iostream>
#include <cassert>
#include <memory>
#include <fstream>


#define all blitz::Range::all()


void onefluid_sim :: run_sim (settingsfile SF)
{
	/*
	 *	Run a simulation in single-fluid mode. In this mode interface tracking and ghost
	 *	fluid methods are turned off. This is a useful test of single-material Riemann
	 *	solvers and Euler solvers.
	 */
	
	std::shared_ptr<eos_base> eos;
	std::shared_ptr<singlefluid_RS_base> RS;
	std::shared_ptr<flow_solver_base> FS;
	
	fluid_state_array statearr (construct_initialise_onefluid(SF, eos, RS, FS));
	fluid_state_array tempstatearr (statearr.copy());

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;

	statearr.output_to_file(SF.basename + std::to_string(numsteps) + ".dat");

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;


	while (t < SF.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		dt = compute_dt(CFL, statearr, SF.T, t);

		FS->single_fluid_update(statearr, tempstatearr, dt);

		statearr.CV = tempstatearr.CV;
		statearr.apply_BCs();
		
		numsteps++;
		t += dt;
		if (SF.output) statearr.output_to_file(SF.basename + std::to_string(numsteps) + ".dat");
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}

	
	
	statearr.output_to_file(SF.basename + "final.dat");
	output_errornorms_to_file(statearr, SF);
	output_cellwise_error(statearr, SF);

	std::cout << "[" << SF.basename << "] Simulation complete." << std::endl;
}




double onefluid_sim :: compute_dt (

	double CFL, 
	fluid_state_array& state, 
	double T, 
	double t
)
{
	/*
	 *	Compute largest stable time step using fluid velocity
	 */

	double maxu = 0.0;

	for (int i=state.array.numGC; i<state.array.length + state.array.numGC; i++)
	{
		maxu = std::max(fabs(state.CV(i,1)/state.CV(i,0)) + state.eos->a(state.CV(i,all)), maxu);
	}

	double dt = CFL*state.array.dx/maxu;

	if (t + dt > T) dt = T - t;

	return dt;
}





void twofluid_sim :: run_sim (settingsfile SF)
{
	/*
	 *	Run a simulation in two fluid mode using the ghost fluid method.
	 */
	
	std::shared_ptr<eos_base> eos1;
	std::shared_ptr<eos_base> eos2;
	std::shared_ptr<singlefluid_RS_base> RS_pure;
	std::shared_ptr<multimat_RS_base> RS_mixed;
	std::shared_ptr<flow_solver_base> FS;
	std::shared_ptr<GFM_base> GFM;
	fluid_state_array statearr1;
	fluid_state_array statearr2;
	levelset_array ls;
	
	construct_initialise_twofluid(SF,eos1,eos2,RS_pure,RS_mixed,FS,GFM,statearr1,statearr2,ls);
	
	fluid_state_array tempstatearr1 (statearr1.copy());
	fluid_state_array tempstatearr2 (statearr2.copy());

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;

	output_endoftimestep(0, SF, statearr1, statearr2, ls);

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;


	while (t < SF.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		dt = compute_dt(CFL, SF.T, t, statearr1, statearr2, ls);
	
		GFM->set_ghost_cells(statearr1, statearr2, ls, RS_mixed);
		FS->single_fluid_update(statearr1, tempstatearr1, dt);
		FS->single_fluid_update(statearr2, tempstatearr2, dt);
		ls.advection_step(dt, GFM->extension_interface_velocity);

		statearr1.CV = tempstatearr1.CV;
		statearr1.apply_BCs();
		statearr2.CV = tempstatearr2.CV;
		statearr2.apply_BCs();
		
		numsteps++;
		t += dt;
		if (SF.output) output_endoftimestep(numsteps, SF, statearr1, statearr2, ls);
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}

	output_endofsimulation(numsteps, SF, statearr1, statearr2, ls);

	std::cout << "[" << SF.basename << "] Simulation complete." << std::endl;
}



double twofluid_sim :: compute_dt (
	
	double CFL, 
	double T, 
	double t, 
	fluid_state_array& state1, 
	fluid_state_array& state2, 
	levelset_array& ls
)
{
	/*
	 *	Compute the largest stable time step in the two fluid sim
	 *	by looking at real velocities
	 */
	
	double maxu = 0.0;

	for (int i=state1.array.numGC; i<state1.array.length + state1.array.numGC; i++)
	{
		double u1 = fabs(state1.CV(i,1)/state1.CV(i,0)) + state1.eos->a(state1.CV(i,all));
		double u2 = fabs(state2.CV(i,1)/state2.CV(i,0)) + state2.eos->a(state2.CV(i,all));
		
		if (ls.linear_interpolation(state1.array.cellcentre_coord(i)) <= 0.0)
		{
			maxu = std::max(maxu, u1);
		}
		else
		{
			maxu = std::max(maxu, u2);
		}
	}

	double dt = CFL*state1.array.dx/maxu;

	if (t + dt > T) dt = T - t;

	return dt;
}


void twofluid_sim :: output_endoftimestep (

	int numsteps, 
	settingsfile& SF, 
	fluid_state_array& state1, 
	fluid_state_array& state2, 
	levelset_array& ls
)
{
	/*
	 *	All outputs to be perfomed upon completion of time step
	 */
	
	state1.output_to_file(SF.basename + "fluid1_" + std::to_string(numsteps) + ".dat");
	state2.output_to_file(SF.basename + "fluid2_" + std::to_string(numsteps) + ".dat");
	ls.output_to_file(SF.basename + "ls_" + std::to_string(numsteps) + ".dat");
	output_realfluidonly(std::to_string(numsteps), SF, state1, state2, ls);
}
	

void twofluid_sim :: output_endofsimulation (

	int numsteps, 
	settingsfile& SF, 
	fluid_state_array& state1, 
	fluid_state_array& state2, 
	levelset_array& ls
)
{
	/*
	 *	All outputs to be performed upon completion of simulation
	 */

	state1.output_to_file(SF.basename + "fluid1_final.dat");
	state2.output_to_file(SF.basename + "fluid2_final.dat");
	ls.output_to_file(SF.basename + "ls_final.dat");
	output_realfluidonly("final", SF, state1, state2, ls);

	fluid_state_array realstate (state1.copy());
	for (int i=state1.array.numGC; i<state1.array.length + state1.array.numGC; i++)
	{
		double x = state1.array.cellcentre_coord(i);

		if (ls.linear_interpolation(x) <= 0.0)
		{
			realstate.CV(i,all) = state1.CV(i,all);
		}
		else
		{
			realstate.CV(i,all) = state2.CV(i,all);
		}
	}
	output_errornorms_to_file(realstate, SF);
	output_cellwise_error(realstate, SF);
}


void twofluid_sim :: output_realfluidonly (

	std::string number, 
	settingsfile& SF, 
	fluid_state_array& state1, 
	fluid_state_array& state2, 
	levelset_array& ls
)
{
	/*
	 *	Output only the real fluid
	 */

	std::ofstream outfile1;
	std::ofstream outfile2;
	outfile1.open(SF.basename + "realfluid1_" + number + ".dat");
	outfile2.open(SF.basename + "realfluid2_" + number + ".dat");

	for (int i=state1.array.numGC; i<state1.array.length + state1.array.numGC; i++)
	{
		double x = state1.array.cellcentre_coord(i);

		if (ls.linear_interpolation(x) <= 0.0)
		{
			outfile1 << x << " " << state1.CV(i,0) << " " << state1.CV(i,1)/state1.CV(i,0) << " " 
				<< state1.CV(i,2) << " " << state1.eos->p(state1.CV(i,all)) << " " << specific_ie_cv(state1.CV(i,all)) << std::endl;
		}
		else
		{
			outfile2 << x << " " << state2.CV(i,0) << " " << state2.CV(i,1)/state2.CV(i,0) << " " 
				<< state2.CV(i,2) << " " << state2.eos->p(state2.CV(i,all)) << " " << specific_ie_cv(state2.CV(i,all)) << std::endl;
		}
	}
}
