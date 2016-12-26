/*
 *	DESCRIPTION:	This file contains a number of tests for various bits of the code.
 *			It can be compiled by running 'make tester'.
 * 
 */



#include "data_storage.hpp"
#include "eos.hpp"
#include "construct_initialise.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "vfield.hpp"
#include "levelset_advection.hpp"
#include "ghost_fluid_method.hpp"
#include <memory>
#include <iostream>







void levelsetadvection_test1();

void ghost_fluid_method_test1();

void full_idealgas_test1();

void test_settingsfile();

void test_M_HLLC_solution ();




int main()
{
	//onefluid_testcase1_godunov_HLLC();
	//onefluid_testcase2_godunov_HLLC();
	//levelsetadvection_test1();
	//ghost_fluid_method_test1();
	//full_idealgas_test1();
	test_M_HLLC_solution();
}








void levelsetadvection_test1 ()
{
	// Simulation parameters

	arrayinfo array;
	array.length = 100;
	array.x0 = 0.0;
	array.dx = 0.01;
	array.numGC = 2;
	array.leftBC = transmissive;
	array.rightBC = transmissive;
	
	ls_IC_type IC = T1;
	double T = 0.5;
	double dt = 0.75*array.dx;

	std::string basename = "./test/testoutput/lsadvection_T1_";


	// Create objects

	std::shared_ptr<vfield_base> vfield = std::make_shared<vfield_test1>();

	levelset_array ls (array);
	levelset_array newls (array);

	std::shared_ptr<levelset_advection_base> lsadvection = std::make_shared<levelset_advection_firstorderup>();


	// Initial conditions

	initialise_levelset (ls, IC);
	ls.apply_BCs();
	int num_iter = 0;
	ls.output_to_file(basename + std::to_string(num_iter) + ".dat");

	std::cout << "[levelsetadvection_test1] Initialisation complete" << std::endl;


	// Time iterations

	double t = 0.0;

	while (t < T)
	{
		lsadvection->advection_operator(ls, newls, vfield, dt);

		ls.phi = newls.phi;
		ls.apply_BCs();

		num_iter++;
		t += dt;
		ls.output_to_file(basename + std::to_string(num_iter) + ".dat");

		std::cout << "[levelsetadvection_test1] Time step " << num_iter << " complete. t = "
				<< t << std::endl;
	}

	std::cout << "[levelsetadvection_test1] Test complete!" << std::endl;
}






void ghost_fluid_method_test1()
{
	// Array parameters

	arrayinfo statearray;
	statearray.length = 100;
	statearray.x0 = 0.0;
	statearray.dx = 0.01;
	statearray.numGC = 4;
	statearray.leftBC = transmissive;
	statearray.rightBC = transmissive;

	arrayinfo lsarray;
	lsarray.length = 500;
	lsarray.x0 = 0.0;
	lsarray.dx = 0.002;
	lsarray.numGC = 2;
	lsarray.leftBC = transmissive;
	lsarray.rightBC = transmissive;

	double gamma1 = 1.4;
	double gamma2 = 1.6;

	ls_IC_type lsIC = T1;
	IC_type IC = TC1;
	
	std::string basename = "./test/testoutput/GFM_original_";


	// Objects used for this test

	std::shared_ptr<ghost_fluid_method_base> GFM = std::make_shared<original_GFM>(lsarray);

	levelset_array ls (lsarray);

	std::shared_ptr<eos_base> eos1 = std::make_shared<eos_idealgas>(gamma1);
	std::shared_ptr<eos_base> eos2 = std::make_shared<eos_idealgas>(gamma2);

	twofluid_array states (statearray,eos1,eos2);

	std::shared_ptr<riemann_solver_base> RS = std::make_shared<HLLC_riemann_solver_idealgas>();


	// Set initial conditions and output

	initialise_levelset(ls,lsIC);
	initialise_twofluid(states,IC);

	ls.output_to_file(basename + "lsinitial.dat");
	states.output_to_file(basename + "statesinitial.dat");


	// Set ghost cells and then output
	
	GFM->set_ghost_cells(states,ls,RS);

	ls.output_to_file(basename + "lsafterGFM.dat");
	states.output_to_file(basename + "statesafterGFM.dat");


	// Set boundary conditions and then output

	ls.apply_BCs();
	states.apply_BCs();

	ls.output_to_file(basename + "lsafterBC.dat");
	states.output_to_file(basename + "statesafterBC.dat");
	
}





void full_idealgas_test1()
{
	// Array parameters

	arrayinfo statearray;
	statearray.length = 100;
	statearray.x0 = 0.0;
	statearray.dx = 0.01;
	statearray.numGC = 1;
	statearray.leftBC = transmissive;
	statearray.rightBC = transmissive;

	arrayinfo lsarray;
	lsarray.length = 100;
	lsarray.x0 = 0.0;
	lsarray.dx = 0.01;
	lsarray.numGC = 2;
	lsarray.leftBC = transmissive;
	lsarray.rightBC = transmissive;

	double gamma1 = 1.4;
	double gamma2 = 1.4;

	ls_IC_type lsIC = T1;
	IC_type IC = TC1;
	
	std::string basename = "./test/testoutput/GFM_isobaricfix_GOD_HLLC_TC1_";
	std::string baselsname = basename + "ls_";

	double T = 0.25;
	double CFL0 = 0.8;


	// Objects used for this test

	std::shared_ptr<ghost_fluid_method_base> GFM = std::make_shared<isobaric_fix_GFM>(lsarray);

	levelset_array ls (lsarray);
	levelset_array newls (lsarray);

	std::shared_ptr<eos_base> eos1 = std::make_shared<eos_idealgas>(gamma1);
	std::shared_ptr<eos_base> eos2 = std::make_shared<eos_idealgas>(gamma2);

	twofluid_array states (statearray,eos1,eos2);
	twofluid_array newstates (statearray,eos1,eos2);

	std::shared_ptr<riemann_solver_base> RS = std::make_shared<HLLC_riemann_solver_idealgas>();
	std::shared_ptr<flow_solver_base> FS = std::make_shared<godunov>(RS);

	std::shared_ptr<vfield_base> vfield = std::make_shared<vfield_starstate>(FS,ls,statearray);
	std::shared_ptr<levelset_advection_base> lsadvection = std::make_shared<levelset_advection_firstorderup>();


	// Set initial conditions and output

	initialise_levelset(ls,lsIC);
	initialise_twofluid(states,IC);

	int num_iter = 0;
	double t = 0.0;
	states.output_to_file(basename + std::to_string(num_iter) + ".dat");
	states.output_realfluid_to_file(basename + "real_" + std::to_string(num_iter) + ".dat", ls);
	states.output_conservation_error(basename + "cv_error.dat", ls, t);
	ls.output_to_file(baselsname + std::to_string(num_iter) + ".dat");

	std::cout << "[GFM_original_GOD_HLLC_TC1] Initialisation complete" << std::endl;


	// Time iterations

	double CFL;

	while (t < T)
	{
		if (num_iter <= 5) CFL = 0.2;
		else CFL = CFL0;
		double dt = compute_dt_serial(CFL, states, ls, T, t);

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

		states.output_to_file(basename + std::to_string(num_iter) + ".dat");
		states.output_realfluid_to_file(basename + "real_" + std::to_string(num_iter) + ".dat", ls);
		states.output_conservation_error(basename + "cv_error.dat", ls, t);
		ls.output_to_file(baselsname + std::to_string(num_iter) + ".dat");

		std::cout << "[GFM_original_GOD_HLLC_TC1] Time step number " << num_iter 
			<< " complete. t = " << t << std::endl;

	}

	std::cout << "[GFM_original_GOD_HLLC_TC1] Test complete!" << std::endl;
}






void test_settingsfile()
{
	
}









void test_M_HLLC_solution()
{
	// Set up left and right fluid states

	 blitz::Array<double,1> leftprimitives (3);
	 leftprimitives(0) = 1.0;
	 leftprimitives(1) = 0.0;
	 leftprimitives(2) = 1.0;

	blitz::Array<double,1> rightprimitives (3);
	rightprimitives(0) = 0.125;
	rightprimitives(1) = 0.0;
	rightprimitives(2) = 0.1;


	// Set up left and right eos

	std::shared_ptr<eos_base> eosL = std::make_shared<eos_idealgas>(1.4);
	std::shared_ptr<eos_base> eosR = std::make_shared<eos_idealgas>(1.4);

        blitz::Array<double,1> leftstate (3);
	leftstate(0) = leftprimitives(0);
        leftstate(1) = leftprimitives(0)*leftprimitives(1);
        leftstate(2) = eosL->E(leftprimitives);
							
        blitz::Array<double,1> rightstate (3);
        rightstate(0) = rightprimitives(0);
        rightstate(1) = rightprimitives(0)*rightprimitives(1);
        rightstate(2) = eosR->E(rightprimitives);
							

	// Call RS for interface state
	
	std::shared_ptr<riemann_solver_base> RS = std::make_shared<M_HLLC_riemann_solver>();
	double p_star, u_star, rho_star_L, rho_star_R;
	RS->solve_rp_forinterfaceboundary(leftstate, rightstate, p_star, u_star, rho_star_L, rho_star_R, eosL, eosR);

	
	// Output

	std::cout << "p_star = " << p_star << std::endl;
	std::cout << "u_star = " << u_star << std::endl;
	std::cout << "rho_star_L = " << rho_star_L << std::endl;
	std::cout << "rho_star_R = " << rho_star_R << std::endl;
}
