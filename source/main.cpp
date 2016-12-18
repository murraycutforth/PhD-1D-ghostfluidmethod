/*
 * 	DESCRIPTION:	Main function for 1D ghost fluid code. After reading settings file, the main function
 *			decides to run a simulation in either single fluid, or full level set - ghost fluid
 *			two fluid mode.
 *
 */


#include "data_storage.hpp"
#include "run_sim.hpp"
#include "input.hpp"
#include <iostream>
#include <cassert>
#include <memory>


int main()
{
	std::cout << "Beginning parallel_GFM code.." << std::endl;

	settingsfile SF;
	SF.read_settings_file();

	std::cout << "Settings file loaded." << std::endl;
	std::cout << "Simulation name is " << SF.basename << std::endl;
	

	std::shared_ptr<sim_base> sim;

	if (SF.sim == onefluid)
	{
		sim = std::make_shared<onefluid_sim>();
	}
	//else if (SF.sim == twofluid)
	//{
	//	sim = std::make_shared<twofluid_sim>();
	//}
	else
	{
		assert(!"Invalid sim type");
	}

	sim->run_sim(SF);

	std::cout << "parallel_GFM code complete!" << std::endl;
	
	return 0;
}
