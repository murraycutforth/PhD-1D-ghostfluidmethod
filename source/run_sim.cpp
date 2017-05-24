/*
 *	DESCRIPTION:	File defines classes which implement various simulation modes. Each definition
 *			of the run_sim member function performs all steps required to solve an
 *			initial-boundary value problem for the Euler equations, and output results
 *			to file.
 *
 *
 */


#include "run_sim.hpp"
#include "misc.hpp"
#include "construct_initialise.hpp"
#include "error.hpp"
#include <iostream>
#include <cassert>
#include <memory>
#include <fstream>


#define all blitz::Range::all()


onefluid_sim :: onefluid_sim ()
:
	eos (),
	RS (),
	FS (),
	FL (3),
	FR (3),
	U0 (3),
	Ut (3)
{}


void onefluid_sim :: run_sim (settingsfile SF)
{
	/*
	 *	Run a simulation in single-fluid mode. In this mode interface tracking and ghost
	 *	fluid methods are turned off. This is a useful test of single-material Riemann
	 *	solvers and Euler solvers.
	 */
	
	fluid_state_array statearr (construct_initialise_onefluid(SF, eos, RS, FS));
	fluid_state_array tempstatearr (statearr.copy());

	compute_total_U_onefluid(statearr, U0);
	compute_total_U_onefluid(statearr, Ut);

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;

	statearr.output_to_file(SF.basename + std::to_string(numsteps) + ".dat");
	output_conservation_errors_to_file(Ut, U0, t, SF);

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;


	while (t < SF.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		
		
		/*
		 * First update all stored variables to next time level
		 */
		 
		dt = compute_dt(CFL, statearr, SF.T, t);
		FS->single_fluid_update(statearr, tempstatearr, dt, FL, FR);
		statearr.CV = tempstatearr.CV;
		
		
		/*
		 * Now deal with output and record-keeping
		 */
		
		numsteps++;
		t += dt;
		if (SF.output) statearr.output_to_file(SF.basename + std::to_string(numsteps) + ".dat");
		update_total_U_onefluid(FL, FR, U0, dt);
		compute_total_U_onefluid(statearr, Ut);
		output_conservation_errors_to_file(Ut, U0, t, SF);
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}

	statearr.output_to_file(SF.basename + "final.dat");
	output_errornorms_to_file(statearr, SF);

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





twofluid_sim :: twofluid_sim ()
:
	eos1 (),
	eos2 (),
	RS_pure	(),
	RS_mixed (),
	FS (),
	GFM (),
	newGFM (),
	statearr1 (),
	statearr2 (),
	ls (),
	FL1 (3),
	FR1 (3),
	FL2 (3),
	FR2 (3),
	U0 (3),
	Ut (3)
{}



void twofluid_sim :: run_sim (settingsfile SF)
{
	/*
	 *	Run a simulation in two fluid mode using the ghost fluid method.
	 */
	
	construct_initialise_twofluid(
	
		SF,
		eos1,
		eos2,
		RS_pure,
		RS_mixed,
		FS,
		GFM,
		newGFM,
		statearr1,
		statearr2,
		ls
	);
	
	fluid_state_array tempstatearr1 (statearr1.copy());
	fluid_state_array tempstatearr2 (statearr2.copy());
	levelset_array prev_ls (ls.copy());

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;

	compute_total_U_twofluid(statearr1, statearr2, ls, U0);
	compute_total_U_twofluid(statearr1, statearr2, ls, Ut);

	output_endoftimestep(0, SF, statearr1, statearr2, ls, t);
	output_conservation_errors_to_file(Ut, U0, t, SF);

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;


	while (t < SF.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		
		
		/*
		 * First update all stored variables to next time level
		 */
		
		dt = compute_dt(CFL, SF.T, t, statearr1, statearr2, ls);
	
		GFM->set_ghost_cells(statearr1, statearr2, ls, prev_ls, RS_mixed, dt);
		
		FS->single_fluid_update(statearr1, tempstatearr1, dt, FL1, FR1);
		FS->single_fluid_update(statearr2, tempstatearr2, dt, FL2, FR2);
		ls.advection_step(dt, GFM->extension_interface_velocity, prev_ls);
		
		newGFM->update_state(statearr1, tempstatearr1, statearr2, tempstatearr2, ls, prev_ls, RS_mixed, dt);

		statearr1.CV = tempstatearr1.CV;
		statearr2.CV = tempstatearr2.CV;
		
		
		/*
		 * Now deal with output and record-keeping
		 */
		
		numsteps++;
		t += dt;
		if (SF.output) output_endoftimestep(numsteps, SF, statearr1, statearr2, ls, t);
		update_total_U_twofluid(FL1, FR1, FL2, FR2, ls, U0, dt, statearr1);
		compute_total_U_twofluid(statearr1, statearr2, ls, Ut);
		output_conservation_errors_to_file(Ut, U0, t, SF);
		
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

	for (int i=0; i<=state1.array.length + state1.array.numGC; i++)
	{
		double u1 = fabs(state1.CV(i,1)/state1.CV(i,0)) + state1.eos->a(state1.CV(i,all));
		double u2 = fabs(state2.CV(i,1)/state2.CV(i,0)) + state2.eos->a(state2.CV(i,all));
		
		maxu = std::max(maxu, u1);
		maxu = std::max(maxu, u2);
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
	levelset_array& ls,
	double t
)
{
	/*
	 *	All outputs to be perfomed upon completion of time step
	 */
	
	static double WM0, AM0;
	 
	if (SF.IC == "NE4")
	{
		// Output wall pressure coefficients at current time
		
		if (numsteps % 10 == 0)
		{
		
			double p0 = 1.0;
			 
			double PL = state1.eos->p(state1.CV(state1.array.numGC,all));
			double PR = state1.eos->p(state1.CV(state1.array.numGC + state1.array.length - 1,all));
			
			double WM = 0.0, AM = 0.0;
			for (int i=state1.array.numGC; i<state1.array.length + state1.array.numGC; i++)
			{
				double x = state1.array.cellcentre_coord(i);
				
				double phiL = ls.linear_interpolation(x - 0.5*state1.array.dx);
				double phiR = ls.linear_interpolation(x + 0.5*state1.array.dx);
				
				if (phiL <= 0.0 && phiR <= 0.0)
				{
					AM += state1.array.dx*state1.CV(i,0);
				}
				else if (phiL > 0.0 && phiR > 0.0)
				{
					WM += state1.array.dx*state2.CV(i,0);
				}
				else
				{
					if (phiL <= 0.0)
					{
						AM += (fabs(phiL))*state1.CV(i,0);
						WM += ((state1.array.dx - fabs(phiL)))*state2.CV(i,0);
					}
					else
					{
						WM += (fabs(phiL))*state2.CV(i,0);
						AM += ((state1.array.dx - fabs(phiL)))*state1.CV(i,0);
					}
				}
			}
			
			std::ofstream outfile1, outfile2, outfile3, outfile4;
			if (t==0)
			{ 
				outfile1.open(SF.basename + "_PL.dat");
				outfile2.open(SF.basename + "_PR.dat");
				outfile3.open(SF.basename + "_watermass.dat");
				outfile4.open(SF.basename + "_airmass.dat");
				
				AM0 = AM;
				WM0 = WM;
			}
				
		
			else
			{
				outfile1.open(SF.basename + "_PL.dat", std::fstream::app);
				outfile2.open(SF.basename + "_PR.dat", std::fstream::app);
				outfile3.open(SF.basename + "_watermass.dat", std::fstream::app);
				outfile4.open(SF.basename + "_airmass.dat", std::fstream::app);
			}
			
			outfile1 << t << " " << (PL - p0)/p0 << std::endl;
			outfile2 << t << " " << (PR - p0)/p0 << std::endl;
			outfile3 << t << " " << (WM - WM0)/WM0 << std::endl;
			outfile4 << t << " " << (AM - AM0)/AM0 << std::endl;
		}
	}
	else
	{
		state1.output_to_file(SF.basename + "fluid1_" + std::to_string(numsteps) + ".dat");
		state2.output_to_file(SF.basename + "fluid2_" + std::to_string(numsteps) + ".dat");
		ls.output_to_file(SF.basename + "ls_" + std::to_string(numsteps) + ".dat");
		output_realfluidonly(std::to_string(numsteps), SF, state1, state2, ls);
	}
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

	if (SF.IC != "NE4")
	{
		output_twofluid_errornorms_to_file(state1, state2, ls, SF);
		output_twofluid_cellwise_error(state1, state2, ls, SF);
	}
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







