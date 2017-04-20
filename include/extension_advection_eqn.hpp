#ifndef EXTENSION_ADVECTION_EQN_H
#define EXTENSION_ADVECTION_EQN_H


#include "data_storage.hpp"
#include <blitz/array.h>


void extension_advection_eqn_1D (
	
	levelset_array& ls,
	fluid_state_array& ghoststatearr1,
	fluid_state_array& ghoststatearr2,
	blitz::Array<double,1> vfield,
	int frozenwidth = 1
);


#endif

