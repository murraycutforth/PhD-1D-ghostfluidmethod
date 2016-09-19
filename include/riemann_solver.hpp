#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER



#include "data_storage.hpp"



#include <blitz/array.h>
#include <memory>




class riemann_solver_base {

	public:

	virtual void solve_rp (	blitz::Array<double,1> Lstate, 
				blitz::Array<double,1> Rstate, 
				blitz::Array<double,1> flux,
				double& S_star,
				std::shared_ptr<eos_base> eosL,
				std::shared_ptr<eos_base> eosR) =0;
};






class HLLC_riemann_solver_idealgas : public riemann_solver_base {

	public:

	void solve_rp (	blitz::Array<double,1> Lstate, 
			blitz::Array<double,1> Rstate, 
			blitz::Array<double,1> flux,
			double& S_star,
			std::shared_ptr<eos_base> eosL,
			std::shared_ptr<eos_base> eosR);


	double q_K (double gamma, double p_star, double p_K);
};



#endif
