#include "data_storage.hpp"
#include "run_sim.hpp"



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

	if (SF.sim == serial_onefluid)
	{
		sim = std::make_shared<serial_onefluid_sim>();
	}
	else if (SF.sim == serial_twofluid)
	{
		sim = std::make_shared<serial_twofluid_sim>();
	}
	else
	{
		assert(!"Invalid sim type");
	}

	sim->run_sim(SF);

	std::cout << "parallel_GFM code complete!" << std::endl;
	
	return 0;
}
