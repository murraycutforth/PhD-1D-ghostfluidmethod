#ifndef FLOW_SOLVER
#define FLOW_SOLVER



#include "data_storage.hpp"
#include "eos.hpp"



#include <blitz/array.h>
#include <memory>
#include <algorithm>




class riemann_solver_base;





blitz::Array<double,1> euler_flux (double rho, double u, double P, double E);

blitz::Array<double,1> euler_flux (blitz::Array<double,1> cv, std::shared_ptr<eos_base> eos);






double limited_slope (double beta, double del_L, double del_R);






class flow_solver_base {

	public:

	blitz::Array<double,1> edge_velocity;	// Stores S_star from solution of riemann problem for interface tracking in non M-RS variants
	
	
	std::shared_ptr<riemann_solver_base> rs;


	flow_solver_base (std::shared_ptr<riemann_solver_base> rs);


	virtual void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt) =0;
};





class godunov : public flow_solver_base {

	public:


	godunov (std::shared_ptr<riemann_solver_base> rs);



	void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt);

};




class MUSCL : public flow_solver_base {

	public:

	MUSCL (std::shared_ptr<riemann_solver_base> rs);


	void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt);

	blitz::Array<double,1> MUSCL_slope (double omega, blitz::Array<double,1> delU_L, blitz::Array<double,1> delU_R);

};



class debugger : public flow_solver_base {

	public:

	debugger (std::shared_ptr<riemann_solver_base> rs);

	void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt);
};



#endif
