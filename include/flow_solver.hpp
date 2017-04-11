#ifndef FLOW_SOLVER
#define FLOW_SOLVER


#include "data_storage.hpp"
#include "eos.hpp"
#include <blitz/array.h>
#include <memory>
#include <algorithm>


class singlefluid_RS_base;


class flow_solver_base {

	public:
	
	std::shared_ptr<singlefluid_RS_base> rs;

	flow_solver_base (std::shared_ptr<singlefluid_RS_base> rs);

	virtual void single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt) =0;
};



class godunov : public flow_solver_base {

	public:

	godunov (std::shared_ptr<singlefluid_RS_base> rs);

	void single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt);
};



class MUSCL : public flow_solver_base {

	public:

	MUSCL (std::shared_ptr<singlefluid_RS_base> rs);

	void single_fluid_update (fluid_state_array& oldstate, fluid_state_array& newstate, double dt);

	blitz::Array<double,1> limited_slope (blitz::Array<double,1> del_L, blitz::Array<double,1> del_R);
};


#endif
