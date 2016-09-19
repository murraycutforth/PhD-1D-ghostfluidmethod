#ifndef FLOW_SOLVER
#define FLOW_SOLVER



#include "data_storage.hpp"



#include <blitz/array.h>
#include <memory>




class riemann_solver_base;





blitz::Array<double,1> euler_flux (double rho, double u, double P, double E);





class flow_solver_base {

	public:

	blitz::Array<double,1> edge_velocity;	// Stores S_star from solution of riemann problem for interface trackin
	
	
	std::shared_ptr<riemann_solver_base> rs;


	flow_solver_base (std::shared_ptr<riemann_solver_base> rs);


	virtual void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt) =0;
};





class godunov : public flow_solver_base {

	public:


	godunov (std::shared_ptr<riemann_solver_base> rs);



	void single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt);

};



#endif
