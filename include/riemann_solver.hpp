#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER



#include "data_storage.hpp"



#include <blitz/array.h>
#include <memory>




class riemann_solver_base {

	public:

	/*
	 * This function is called inside the flow solver to obtain fluxes on the (x/t)=0 characteristic
	 * Only called within single fluid update.
	 */

	virtual void solve_rp (	blitz::Array<double,1> Lstate, 
				blitz::Array<double,1> Rstate, 
				blitz::Array<double,1> flux,
				double& S_star,
				std::shared_ptr<eos_base> eos) =0;
	


	/*
	 * This function is called by the ghost fluid method to solve the multimaterial riemann problem.
	 * It finds the interface states (pressure, particle speed, density to L and R)
	 */
	virtual void solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
							blitz::Array<double,1> Rstate,
							double& p_star,
							double& u_star,
							double& rho_star_L,
							double& rho_star_R,
							std::shared_ptr<eos_base> eosL,
							std::shared_ptr<eos_base> eosR ) =0;
};






class HLLC_riemann_solver_idealgas : public riemann_solver_base {

	public:

	void solve_rp (	blitz::Array<double,1> Lstate, 
			blitz::Array<double,1> Rstate, 
			blitz::Array<double,1> flux,
			double& S_star,
			std::shared_ptr<eos_base> eos);

	void solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR );


	double q_K (double gamma, double p_star, double p_K);
};



class M_HLLC_riemann_solver : public riemann_solver_base {

	public:

	void solve_rp (	blitz::Array<double,1> Lstate,
			blitz::Array<double,1> Rstate,
			blitz::Array<double,1> flux,
			double& S_star,
			std::shared_ptr<eos_base> eos);
	
	void solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR );

	double mu (double fL, double fR, double rhoL, double rhoR);

};



class exact_riemann_solver_idealgas : public riemann_solver_base {

	public:

	void solve_rp (	blitz::Array<double,1> Lstate,
			blitz::Array<double,1> Rstate,
			blitz::Array<double,1> flux,
			double& S_star,
			std::shared_ptr<eos_base> eos);
	
	void solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR );
};




#endif
