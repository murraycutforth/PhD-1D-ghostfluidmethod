#ifndef CONSTRUCT_INITIALISE_H
#define CONSTRUCT_INITIALISE_H


#include "input.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "data_storage.hpp"
#include <memory>
#include <blitz/array.h>


fluid_state_array construct_initialise_onefluid ( 

	settingsfile& SF, 
	std::shared_ptr<eos_base>& eos, 
	std::shared_ptr<singlefluid_RS_base>& RS, 
	std::shared_ptr<flow_solver_base>& FS
);


void set_piecewiseconstant_ICs (

	blitz::Array<double,1> leftprimitives,
	blitz::Array<double,1> rightprimitives,
	double discontinuitylocation,
	fluid_state_array& statearr
);


#endif
