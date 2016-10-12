#include "riemann_solver.hpp"
#include "flow_solver.hpp"
#include "eos.hpp"
#include "exact_RS_idealgas.hpp"



#include <cassert>
#include <cmath>
#include <algorithm>






void HLLC_riemann_solver_idealgas ::  solve_rp (	blitz::Array<double,1> Lstate, 
							blitz::Array<double,1> Rstate, 
							blitz::Array<double,1> flux,
							double& S_star,
							std::shared_ptr<eos_base> eos)
{
	assert(Lstate.extent(blitz::firstDim) == 3);
	assert(Rstate.extent(blitz::firstDim) == 3);
	assert(flux.extent(blitz::firstDim) == 3);


	// Recover the primitive variables

	double rho_L = Lstate(0);
	double u_L = Lstate(1)/Lstate(0);
	double p_L = eos->p(Lstate);

	double rho_R = Rstate(0);
	double u_R = Rstate(1)/Rstate(0);
	double p_R = eos->p(Rstate);


	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));
	assert(eos->get_eos_type() == ideal);
	assert(eos->get_gamma() == eos->get_gamma());



	// Estimate pressure using PVRS scheme (linearising around arithmetic average state)

	double a_L = eos->a(Lstate);
	double a_R = eos->a(Rstate);
	double p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R-u_L)*0.5*(rho_L+rho_R)*0.5*(a_L+a_R);
	double p_star = std::max(0.0,p_pvrs);


	// Pressure-based wave speed estimates

	double S_L = u_L - a_L*q_K(eos->get_gamma(), p_star, p_L);
	double S_R = u_R + a_R*q_K(eos->get_gamma(), p_star, p_R);
	S_star = (p_R-p_L+rho_L*u_L*(S_L-u_L)-rho_R*u_R*(S_R-u_R))/(rho_L*(S_L-u_L)-rho_R*(S_R-u_R));


	// Return appropriate HLLC flux

	if (0.0 <= S_L)
	{
		flux = euler_flux(rho_L, u_L, p_L, Lstate(2));
	}
	else if (S_L <= 0.0 && 0.0 <= S_star)
	{
		double factor_L = rho_L*((S_L - u_L)/(S_L - S_star));

		blitz::Array<double,1> Lstar_state (3);
		Lstar_state(0) = factor_L;
		Lstar_state(1) = factor_L*S_star;
		Lstar_state(2) = factor_L*((Lstate(2)/rho_L) + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))));

		flux = euler_flux(rho_L, u_L, p_L, Lstate(2)) + S_L*(Lstar_state - Lstate);
	}
	else if (S_star <= 0.0 && 0.0 <= S_R)
	{
		double factor_R = rho_R*((S_R - u_R)/(S_R - S_star));

		blitz::Array<double,1> Rstar_state (3);
		Rstar_state(0) = factor_R;
		Rstar_state(1) = factor_R*S_star;
		Rstar_state(2) = factor_R*((Rstate(2)/rho_R) + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))));

		flux = euler_flux(rho_R, u_R, p_R, Rstate(2)) + S_R*(Rstar_state - Rstate);
	}
	else
	{
		assert(S_R <= 0.0);
		flux = euler_flux(rho_R, u_R, p_R, Rstate(2));
	}


}





double HLLC_riemann_solver_idealgas :: q_K (double gamma, double p_star, double p_K)
{
	/*
	 * This function is used in the pressure-based estimate of the wave speeds.
	 * It differentiates between shock and rarefaction waves. See Toro p330.
	 */

	if (p_star <= p_K)
	{
		return 1.0;
	}
	else
	{
		double rvs = 1.0 + ((gamma+1.0)/(2.0*gamma))*((p_star/p_K) - 1.0);
		assert(rvs >= 0.0);
		return sqrt(rvs);
	}
}



void HLLC_riemann_solver_idealgas :: solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR )
{
	// This RS is not a mixed riemann solver so this method should not be called

	assert(false);
}















void M_HLLC_riemann_solver ::  solve_rp (	blitz::Array<double,1> Lstate, 
						blitz::Array<double,1> Rstate, 
						blitz::Array<double,1> flux,
						double& S_star,
						std::shared_ptr<eos_base> eos)
{
	assert(Lstate.extent(blitz::firstDim) == 3);
	assert(Rstate.extent(blitz::firstDim) == 3);
	assert(flux.extent(blitz::firstDim) == 3);


	// Recover the primitive variables

	double rho_L = Lstate(0);
	double u_L = Lstate(1)/Lstate(0);
	double p_L = eos->p(Lstate);
	double a_L = eos->a(Lstate);

	double rho_R = Rstate(0);
	double u_R = Rstate(1)/Rstate(0);
	double p_R = eos->p(Rstate);
	double a_R = eos->a(Rstate);


	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));



	// Estimate wave speeds using M-HLLC method [Hu et al. 2009]

	double Tau_L = eos->get_Tau(Lstate);
	double Psi_L = eos->get_Psi(Lstate);

	double Tau_R = eos->get_Tau(Rstate);
	double Psi_R = eos->get_Psi(Rstate);


	// Calculate the Roe-averages

	double u_tilde = mu(u_L, u_R, rho_L, rho_R);
	double porho_tilde = mu(p_L/rho_L, p_R/rho_R, rho_L, rho_R) + 0.5*((u_R-u_L)/(sqrt(rho_L)+sqrt(rho_R)))*((u_R-u_L)/(sqrt(rho_L)+sqrt(rho_R)));
	double Tau_tilde = mu(Tau_L, Tau_R, rho_L, rho_R);
	double Psi_tilde = mu(Psi_L, Psi_R, rho_L, rho_R);
	double a_tilde_sq = Psi_tilde + Tau_tilde*porho_tilde;

	assert(a_tilde_sq >= 0.0);
	double a_tilde = sqrt(a_tilde_sq);


	// Estimate wave speeds

	double S_L = std::min(u_L - a_L, u_tilde - a_tilde);
	double S_R = std::max(u_R + a_R, u_tilde + a_tilde);
	


	S_star = (p_R-p_L+rho_L*u_L*(S_L-u_L)-rho_R*u_R*(S_R-u_R))/(rho_L*(S_L-u_L)-rho_R*(S_R-u_R));


	// Return appropriate HLLC flux

	if (0.0 <= S_L)
	{
		flux = euler_flux(rho_L, u_L, p_L, Lstate(2));
	}
	else if (S_L <= 0.0 && 0.0 <= S_star)
	{
		double factor_L = rho_L*((S_L - u_L)/(S_L - S_star));

		blitz::Array<double,1> Lstar_state (3);
		Lstar_state(0) = factor_L;
		Lstar_state(1) = factor_L*S_star;
		Lstar_state(2) = factor_L*((Lstate(2)/rho_L) + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))));

		flux = euler_flux(rho_L, u_L, p_L, Lstate(2)) + S_L*(Lstar_state - Lstate);
	}
	else if (S_star <= 0.0 && 0.0 <= S_R)
	{
		double factor_R = rho_R*((S_R - u_R)/(S_R - S_star));

		blitz::Array<double,1> Rstar_state (3);
		Rstar_state(0) = factor_R;
		Rstar_state(1) = factor_R*S_star;
		Rstar_state(2) = factor_R*((Rstate(2)/rho_R) + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))));

		flux = euler_flux(rho_R, u_R, p_R, Rstate(2)) + S_R*(Rstar_state - Rstate);
	}
	else
	{
		assert(S_R <= 0.0);
		flux = euler_flux(rho_R, u_R, p_R, Rstate(2));
	}


}
	
	
	
	
	
double M_HLLC_riemann_solver :: mu (double fL, double fR, double rhoL, double rhoR)
{
	return (sqrt(rhoL)*fL + sqrt(rhoR)*fR)/(sqrt(rhoL)+sqrt(rhoR));
}



void M_HLLC_riemann_solver :: solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR )
{
	assert(Lstate.extent(blitz::firstDim) == 3);
	assert(Rstate.extent(blitz::firstDim) == 3);


	// Recover the primitive variables

	double rho_L = Lstate(0);
	double u_L = Lstate(1)/Lstate(0);
	double p_L = eosL->p(Lstate);
	double a_L = eosL->a(Lstate);

	double rho_R = Rstate(0);
	double u_R = Rstate(1)/Rstate(0);
	double p_R = eosR->p(Rstate);
	double a_R = eosR->a(Rstate);


	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));



	// Estimate wave speeds using M-HLLC method [Hu et al. 2009]

	double Tau_L = eosL->get_Tau(Lstate);
	double Psi_L = eosL->get_Psi(Lstate);

	double Tau_R = eosR->get_Tau(Rstate);
	double Psi_R = eosR->get_Psi(Rstate);


	// Calculate the Roe-averages

	double u_tilde = mu(u_L, u_R, rho_L, rho_R);
	double porho_tilde = mu(p_L/rho_L, p_R/rho_R, rho_L, rho_R) + 0.5*((u_R-u_L)/(sqrt(rho_L)+sqrt(rho_R)))*((u_R-u_L)/(sqrt(rho_L)+sqrt(rho_R)));
	double Tau_tilde = mu(Tau_L, Tau_R, rho_L, rho_R);
	double Psi_tilde = mu(Psi_L, Psi_R, rho_L, rho_R);
	double a_tilde_sq = Psi_tilde + Tau_tilde*porho_tilde;

	assert(a_tilde_sq >= 0.0);
	double a_tilde = sqrt(a_tilde_sq);


	// Estimate wave speeds

	double S_L = std::min(u_L - a_L, u_tilde - a_tilde);
	double S_R = std::max(u_R + a_R, u_tilde + a_tilde);
	
	
	// Finally calculate interface states

	u_star = (p_R-p_L+rho_L*u_L*(S_L-u_L)-rho_R*u_R*(S_R-u_R))/(rho_L*(S_L-u_L)-rho_R*(S_R-u_R));

	p_star = p_L + rho_L*(u_L - S_L)*(u_L - u_star);
	

	// The following gives unsatisfactory results!
	//rho_star_L = rho_L*(u_L - S_L)/(u_star - S_L);
	//
	//rho_star_R = rho_R*(u_R - S_R)/(u_star - S_R);


	// Instead try finding exact solution for density using approximate u_star and p_star

	if (p_star > p_R)
	{
		// Right shock

		rho_star_R = eosR->postshock_density(p_star, p_R, rho_R);
	}
	else
	{
		// Right rarefaction wave

		rho_star_R = eosR->postrarefaction_density(p_star, p_R, rho_R);
	}
	
	if (p_star > p_L)
	{
		// Left shock

		rho_star_L = eosL->postshock_density(p_star, p_L, rho_L);
	}
	else
	{
		// Left rarefaction wave

		rho_star_L = eosL->postrarefaction_density(p_star, p_L, rho_L);
	}
}















void exact_riemann_solver_idealgas :: solve_rp (	blitz::Array<double,1> Lstate,
			blitz::Array<double,1> Rstate,
			blitz::Array<double,1> flux,
			double& S_star,
			std::shared_ptr<eos_base> eos)
{
	assert(eos->get_eos_type() == ideal);

	blitz::Array<double,1> Lprimitives (3);
	blitz::Array<double,1> Rprimitives (3);

	Lprimitives(0) = Lstate(0);
	Lprimitives(1) = Lstate(1)/Lstate(0);
	Lprimitives(2) = eos->p(Lstate);

	Rprimitives(0) = Rstate(0);
	Rprimitives(1) = Rstate(1)/Rstate(0);
	Rprimitives(2) = eos->p(Rstate);

	exact_rs_idealgas RS (eos->get_gamma(), eos->get_gamma());
	RS.solve_RP(Lprimitives,Rprimitives);

	blitz::Array<double,1> soln (3);
	soln = RS.sample_solution(0.0);
	double E = eos->E(soln);

	S_star = RS.S_STAR;
	flux = euler_flux(soln(0), soln(1), soln(2), E);
}
	





void exact_riemann_solver_idealgas :: solve_rp_forinterfaceboundary (	blitz::Array<double,1> Lstate,
						blitz::Array<double,1> Rstate,
						double& p_star,
						double& u_star,
						double& rho_star_L,
						double& rho_star_R,
						std::shared_ptr<eos_base> eosL,
						std::shared_ptr<eos_base> eosR )
{
	assert(eosL->get_eos_type() == ideal);
	assert(eosR->get_eos_type() == ideal);

	blitz::Array<double,1> Lprimitives (3);
	blitz::Array<double,1> Rprimitives (3);

	Lprimitives(0) = Lstate(0);
	Lprimitives(1) = Lstate(1)/Lstate(0);
	Lprimitives(2) = eosL->p(Lstate);

	Rprimitives(0) = Rstate(0);
	Rprimitives(1) = Rstate(1)/Rstate(0);
	Rprimitives(2) = eosR->p(Rstate);

	exact_rs_idealgas RS (eosL->get_gamma(), eosR->get_gamma());
	RS.solve_RP(Lprimitives,Rprimitives);
	
	p_star = RS.P_STAR;
	u_star = RS.S_STAR;
	rho_star_L = RS.W_STAR_L(0);
	rho_star_R = RS.W_STAR_R(0);

	std::cout << "p_star = " << p_star << std::endl;
	std::cout << "u_star = " << u_star << std::endl;
	std::cout << "rho_star_L = " << rho_star_L << std::endl;
	std::cout << "rho_star_R = " << rho_star_R << std::endl;
}
