#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER


#include "data_storage.hpp"
#include <blitz/array.h>
#include <memory>


class singlefluid_RS_base {

	public:

	virtual void solve_rp (	
	
		blitz::Array<double,1> Lstate, 
		blitz::Array<double,1> Rstate, 
		blitz::Array<double,1> flux,			
		std::shared_ptr<eos_base> eos
	) =0;
};


class multimat_RS_base {

	public:

	virtual void solve_rp_forinterfaceboundary (	
	
		blitz::Array<double,1> Lstate,
		blitz::Array<double,1> Rstate,
		double& p_star,
		double& u_star,
		double& rho_star_L,
		double& rho_star_R,
		std::shared_ptr<eos_base> eosL,
		std::shared_ptr<eos_base> eosR 
	) =0;
};






class HLLC_RS_idealgas : public singlefluid_RS_base {

	public:

	void solve_rp (	
		
		blitz::Array<double,1> Lstate, 
		blitz::Array<double,1> Rstate, 
		blitz::Array<double,1> flux,
		std::shared_ptr<eos_base> eos
	);

	double q_K (double gamma, double p_star, double p_K);
};


class exact_RS_idealgas : public singlefluid_RS_base {

	public:

	void solve_rp (	
	
		blitz::Array<double,1> Lstate, 
		blitz::Array<double,1> Rstate, 
		blitz::Array<double,1> flux,
		std::shared_ptr<eos_base> eos
	);
};






class M_HLLC_RS : public multimat_RS_base {

	public:

	void solve_rp_forinterfaceboundary (	
	
		blitz::Array<double,1> Lstate,
		blitz::Array<double,1> Rstate,
		double& p_star,
		double& u_star,
		double& rho_star_L,
		double& rho_star_R,
		std::shared_ptr<eos_base> eosL,
		std::shared_ptr<eos_base> eosR 
	);

	double mu (double fL, double fR, double rhoL, double rhoR);

};



class exact_RS_multi_idealgas : public multimat_RS_base {

	public:
	
	void solve_rp_forinterfaceboundary (	
	
		blitz::Array<double,1> Lstate,
		blitz::Array<double,1> Rstate,
		double& p_star,
		double& u_star,
		double& rho_star_L,
		double& rho_star_R,
		std::shared_ptr<eos_base> eosL,
		std::shared_ptr<eos_base> eosR 
	);
};




#endif
