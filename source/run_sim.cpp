#include "run_sim.hpp"
#include "eos.hpp"
#include "initial_conditions.hpp"
#include "timestep.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "vfield.hpp"
#include "levelset_advection.hpp"
#include "ghost_fluid_method.hpp"






#include <iostream>
#include <cassert>
#include <memory>











void serial_onefluid_sim :: run_sim (settingsfile SF)
{
	
	arrayinfo statearray;
	statearray.length = SF.length;
	statearray.x0 = SF.x0;
	statearray.dx = SF.dx;
	statearray.numGC = SF.numGC;
	statearray.leftBC = SF.BC_L;
	statearray.rightBC = SF.BC_R;
	
	
	std::shared_ptr<eos_base> eos;
	if (SF.eos1 == ideal) eos = std::make_shared<eos_idealgas>(SF.fluid1_gamma);

	onefluid_array state (statearray, eos);
	onefluid_array newstate (statearray, eos);

	std::shared_ptr<riemann_solver_base> RS;
	if (SF.RS == HLLC) RS = std::make_shared<HLLC_riemann_solver_idealgas>();
	else if (SF.RS == M_HLLC) RS = std::make_shared<M_HLLC_riemann_solver>();

	std::shared_ptr<flow_solver_base> FS;
	if (SF.FS == Godunov) FS = std::make_shared<godunov>(RS);
	else if (SF.FS == MUSCL_FS) FS = std::make_shared<MUSCL>(RS);
	else assert(!"Invalid flow solver IC");

	initialise_onefluid(state, SF.IC);
	state.apply_BCs();

	int num_iter = 0;
	double t = 0.0;
	state.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
	state.output_conservation_error(SF.basename + "cv_error.dat", t);
	
	double CFL;

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.."  << std::endl;
	while (t < SF.T)
	{
		if (num_iter <= 5) CFL = 0.2;
		else CFL = SF.CFL;
		double dt = compute_dt_serial_onefluid(CFL, state, SF.T, t);

		FS->single_fluid_update(state, newstate, dt);


		state.fluid = newstate.fluid;
		state.apply_BCs();


		if (SF.IC == TC1 && SF.T==0.25)
		{
			state.initial_total_cv(1) += (dt)*0.9;
		}

		num_iter++;
		t += dt;
		
		state.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
		state.output_conservation_error(SF.basename + "cv_error.dat", t);

		std::cout << "[" << SF.basename << "] Time step " << num_iter << " complete. t = " << t << std::endl;
	}

}










void serial_twofluid_sim :: run_sim (settingsfile SF)
{
	arrayinfo statearray;
	statearray.length = SF.length;
	statearray.x0 = SF.x0;
	statearray.dx = SF.dx;
	statearray.numGC = SF.numGC;
	statearray.leftBC = SF.BC_L;
	statearray.rightBC = SF.BC_R;

	arrayinfo lsarray;
	lsarray.length = SF.lslength;
	lsarray.x0 = SF.x0;
	lsarray.dx = SF.lsdx;
	lsarray.numGC = SF.lsnumGC;
	lsarray.leftBC = SF.BC_L;
	lsarray.rightBC = SF.BC_R;




	std::shared_ptr<eos_base> eos1;
	if (SF.eos1 == ideal) eos1 = std::make_shared<eos_idealgas>(SF.fluid1_gamma);

	std::shared_ptr<eos_base> eos2;
	if (SF.eos2 == ideal) eos2 = std::make_shared<eos_idealgas>(SF.fluid2_gamma);

	twofluid_array states (statearray, eos1, eos2);
	twofluid_array newstates (statearray, eos1, eos2);

	levelset_array ls (lsarray);
	levelset_array newls (lsarray);
	
	std::shared_ptr<ghost_fluid_method_base> GFM;
	if (SF.GFM == Original) GFM = std::make_shared<original_GFM>(statearray);
	else if (SF.GFM == Isobaricfix) GFM = std::make_shared<isobaric_fix_GFM>(statearray);
	else if (SF.GFM == Real) GFM = std::make_shared<rGFM>(statearray);
	
	std::shared_ptr<riemann_solver_base> RS;
	if (SF.RS == HLLC) RS = std::make_shared<HLLC_riemann_solver_idealgas>();
	else if (SF.RS == M_HLLC) RS = std::make_shared<M_HLLC_riemann_solver>();

	std::shared_ptr<flow_solver_base> FS;
	if (SF.FS == Godunov) FS = std::make_shared<godunov>(RS);
	else if (SF.FS == MUSCL_FS) FS = std::make_shared<MUSCL>(RS);
	

	// NOt part of settings file yet
	std::shared_ptr<vfield_base> vfield = std::make_shared<vfield_starstate>(FS, ls, statearray);
	std::shared_ptr<levelset_advection_base> lsadvection = std::make_shared<levelset_advection_firstorderup>();


	initialise_levelset(ls, SF.lsIC);
	initialise_twofluid(states, SF.IC);

	int num_iter = 0;
	double t = 0.0;

	states.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
	states.output_realfluid_to_file(SF.basename + std::to_string(num_iter) + ".dat", ls);
	states.output_conservation_error(SF.basename + "cv_error.dat", ls, t);
	ls.output_to_file(SF.basename + "ls_" + std::to_string(num_iter) + ".dat");

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;

	double CFL;

	while (t < SF.T)
	{
		if (num_iter <= 5) CFL = 0.2;
		else CFL = SF.CFL;
		double dt = compute_dt_serial(CFL, states, ls, SF.T, t);

		advance_time_level (	dt,
					states,
					newstates,
					ls,
					newls,
					GFM,
					FS,
					FS,
					lsadvection,
					vfield,
					RS);

		newstates.apply_BCs();
		newls.apply_BCs();

		t += dt;
		num_iter++;

		states.fluid1 = newstates.fluid1;
		states.fluid2 = newstates.fluid2;
		ls.phi = newls.phi;

		
		// Manually update the total conservation stuff

		if (SF.IC == TC1 && SF.T==0.25)
		{
			states.initial_total_cv(1) += dt*0.9;
		}

		states.output_to_file(SF.basename + std::to_string(num_iter) + ".dat");
		states.output_realfluid_to_file(SF.basename + std::to_string(num_iter) + ".dat", ls);
		states.output_conservation_error(SF.basename + "cv_error.dat", ls, t);
		ls.output_to_file(SF.basename + "ls_" + std::to_string(num_iter) + ".dat");

		std::cout << "[" << SF.basename << "] Time step number " << num_iter 
			<< " complete. t = " << t << std::endl;
	}

std::cout << "[" << SF.basename << "] Simulation complete!" << std::endl;
}
