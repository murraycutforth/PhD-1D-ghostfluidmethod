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
#include "construct_initialise.hpp"
#include "error.hpp"
#include <iostream>
#include <cassert>
#include <memory>




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
		maxu = std::max(fabs(state.CV(i,1)/state.CV(i,0)) + state.eos->a(state.CV(i,blitz::Range::all())), maxu);
	}

	double dt = CFL*state.array.dx/maxu;

	if (t + dt > T) dt = T - t;

	return dt;
}





// old code below!
//
//
//
//
//
//
//void serial_twofluid_sim :: run_sim (settingsfile SF)
//{
//	arrayinfo statearray;
//	statearray.length = SF.length;
//	statearray.x0 = SF.x0;
//	statearray.dx = SF.dx;
//	statearray.numGC = SF.numGC;
//	statearray.leftBC = SF.BC_L;
//	statearray.rightBC = SF.BC_R;
//
//	arrayinfo lsarray;
//	lsarray.length = SF.lslength;
//	lsarray.x0 = SF.x0;
//	lsarray.dx = SF.lsdx;
//	lsarray.numGC = SF.lsnumGC;
//	lsarray.leftBC = SF.BC_L;
//	lsarray.rightBC = SF.BC_R;
//
//
//
//
//	std::shared_ptr<eos_base> eos1;
//	if (SF.eos1 == ideal) eos1 = std::make_shared<eos_idealgas>(SF.fluid1_gamma);
//
//	std::shared_ptr<eos_base> eos2;
//	if (SF.eos2 == ideal) eos2 = std::make_shared<eos_idealgas>(SF.fluid2_gamma);
//
//	twofluid_array states (statearray, eos1, eos2);
//	twofluid_array newstates (statearray, eos1, eos2);
//
//	levelset_array ls (lsarray);
//	levelset_array newls (lsarray);
//	
//	std::shared_ptr<ghost_fluid_method_base> GFM;
//	if (SF.GFM == Original) GFM = std::make_shared<original_GFM>(lsarray);
//	else if (SF.GFM == Isobaricfix) GFM = std::make_shared<isobaric_fix_GFM>(lsarray);
//	else if (SF.GFM == Real) GFM = std::make_shared<rGFM>(lsarray);
//	
//	std::shared_ptr<riemann_solver_base> RS;
//	if (SF.RS == HLLC) RS = std::make_shared<HLLC_riemann_solver_idealgas>();
//	else if (SF.RS == M_HLLC) RS = std::make_shared<M_HLLC_riemann_solver>();
//	else if (SF.RS == exact_idealgas) RS = std::make_shared<exact_riemann_solver_idealgas>();
//
//	std::shared_ptr<flow_solver_base> FS;
//	if (SF.FS == Godunov) FS = std::make_shared<godunov>(RS);
//	else if (SF.FS == MUSCL_FS) FS = std::make_shared<MUSCL>(RS);
//	
//
//	// NOt part of settings file yet
//	std::shared_ptr<vfield_base> vfield;
//	if (SF.GFM == Real)
//	{
//		vfield = std::make_shared<vfield_mixedRPsolution>(GFM, ls);
//	}
//	else
//	{
//		vfield = std::make_shared<vfield_starstate>(FS, ls, statearray);
//	}
//	std::shared_ptr<levelset_advection_base> lsadvection = std::make_shared<levelset_advection_firstorderup>();
//
//
//	initialise_levelset(ls, SF.lsIC);
//	initialise_twofluid(states, SF.IC);
//
//	int num_iter = 0;
//	double t = 0.0;
//
//	states.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
//	states.output_realfluid_to_file(SF.basename + std::to_string(num_iter) + ".dat", ls);
//	states.output_conservation_error(SF.basename + "cv_error.dat", ls, t);
//	ls.output_to_file(SF.basename + "ls_" + std::to_string(num_iter) + ".dat");
//
//	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;
//
//	double CFL;
//
//	while (t < SF.T)
//	{
//		if (num_iter <= 5) CFL = 0.2;
//		else CFL = SF.CFL;
//		double dt = compute_dt_serial(CFL, states, ls, SF.T, t);
//
//		advance_time_level (	dt,
//					states,
//					newstates,
//					ls,
//					newls,
//					GFM,
//					FS,
//					FS,
//					lsadvection,
//					vfield,
//					RS);
//
//		newstates.apply_BCs();
//		newls.apply_BCs();
//
//		t += dt;
//		num_iter++;
//
//		states.fluid1 = newstates.fluid1;
//		states.fluid2 = newstates.fluid2;
//		ls.phi = newls.phi;
//
//		
//		// Manually update the total conservation stuff
//
//		if (SF.IC == TC1 && SF.T==0.25)
//		{
//			states.initial_total_cv(1) += dt*0.9;
//		}
//
//		states.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
//		states.output_realfluid_to_file(SF.basename + std::to_string(num_iter) + ".dat", ls);
//		states.output_conservation_error(SF.basename + "cv_error.dat", ls, t);
//		ls.output_to_file(SF.basename + "ls_" + std::to_string(num_iter) + ".dat");
//
//		std::cout << "[" << SF.basename << "] Time step number " << num_iter 
//			<< " complete. t = " << t << std::endl;
//	}
//
//std::cout << "[" << SF.basename << "] Simulation complete!" << std::endl;
//}
