#ifndef FLOW_SOLVER
#define FLOW_SOLVER


#include "data_storage.hpp"
#include "eos.hpp"
#include <blitz/array.h>
#include <memory>
#include <algorithm>


class singlefluid_RS_base;



blitz::Array<double,1> euler_flux (double rho, double u, double P, double E);

blitz::Array<double,1> euler_flux (blitz::Array<double,1> cv, std::shared_ptr<eos_base> eos);



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

	blitz::Array<double,1> MUSCL_slope (double omega, blitz::Array<double,1> delU_L, blitz::Array<double,1> delU_R);

	double limited_slope (double beta, double del_L, double del_R);
};


#endif
