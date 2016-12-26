#ifndef ERROR_H
#define ERROR_H


#include "input.hpp"
#include "data_storage.hpp"
#include <blitz/array.h>


blitz::Array<double,2> get_cellwise_error (
	
	fluid_state_array& fluid1,
	settingsfile& SF
);


void get_density_errornorms (

	blitz::Array<double,2> cellwise_error,
	double& L1error,
	double& Linferror
);


void output_errornorms_to_file (

	fluid_state_array& fluid1,
	settingsfile& SF
);
	

#endif
