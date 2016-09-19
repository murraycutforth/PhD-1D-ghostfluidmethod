#ifndef TIMESTEP
#define TIMESTEP



#include "data_storage.hpp"




#include <memory>
#include <cmath>




class ghost_fluid_method_base;
class flow_solver_base;
class levelset_advection_base;
class vfield_base;
class riemann_solver_base;



double compute_dt_serial (double CFL, twofluid_array& states, levelset_array& ls, double T, double t);



double compute_dt_serial_onefluid (double CFL, onefluid_array& state, double T, double t);



void advance_time_level (	double dt, 
				twofluid_array& states, 
				twofluid_array& newstates, 
				levelset_array& ls,
				levelset_array& newls,
				std::shared_ptr<ghost_fluid_method_base> gfm,
				std::shared_ptr<flow_solver_base> fs1,
				std::shared_ptr<flow_solver_base> fs2,
				std::shared_ptr<levelset_advection_base> lsadvection,
				std::shared_ptr<vfield_base> vfield,
				std::shared_ptr<riemann_solver_base> RS);




#endif
