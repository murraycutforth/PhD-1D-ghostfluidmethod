#ifndef CONSTRUCT_INITIALISE_H
#define CONSTRUCT_INITIALISE_H


#include "input.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "ghost_fluid_method.hpp"
#include "new_ghost_fluid_method.hpp"
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


void construct_initialise_twofluid (
	
	settingsfile& SF, 
	std::shared_ptr<eos_base>& eos1,
	std::shared_ptr<eos_base>& eos2, 
	std::shared_ptr<singlefluid_RS_base>& RS_pure,
	std::shared_ptr<multimat_RS_base>& RS_mixed, 
	std::shared_ptr<flow_solver_base>& FS,
	std::shared_ptr<GFM_base>& GFM,
	std::shared_ptr<newGFM_base>& newGFM,
	fluid_state_array& statearr1,
	fluid_state_array& statearr2,
	levelset_array& ls
);


void set_singlediscontinuity_ls_IC (

	double discontinuitylocation, 
	int parity, 
	levelset_array& ls
);

#endif
