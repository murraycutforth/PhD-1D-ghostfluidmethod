/*
 *	DESCRIPTION: 	This file contains implementations of various Riemann solvers for the Euler equation. The Riemann 
 *			problem is the IVP of two piecewise-constant states separated by a discontinuity at x=0. The single
 *			material Riemann solvers are called by flow solvers to compute the intercell flux in the update
 *			of a single material array. The multi-material Riemann solvers are called by ghost fluid methods
 *			in order to compute the pressure, densities, and velocity across the interface.
 *
 *	CITATIONS: 	E Toro - "Riemann solvers and numerical methods for fluid dynamics: A practical introduction" - 1999
 *			X Hu, N Adams - "On the HLLC Riemann solver for interface interaction in compressible multi-fluid flow" - 2009
 *
 */


#include "riemann_solver.hpp"
#include "flow_solver.hpp"
#include "eos.hpp"
#include "exact_RS_idealgas.hpp"
#include "misc.hpp"
#include <cassert>
#include <string>
#include <cmath>
#include <algorithm>



// ========================================= SINGLE MATERIAL RIEMANN SOLVERS ================================== //


void HLLC_RS_idealgas ::  solve_rp (	

	blitz::Array<double,1> Lstate, 
	blitz::Array<double,1> Rstate, 
	blitz::Array<double,1> flux,
	std::shared_ptr<eos_base> eos
)
{
	/*
	 *	The HLLC approximate Riemann solver for ideal gas. Flux across x/t=0
	 *	characteristic returned.
	 */

	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));
	assert(eos->get_eos_type() == "ideal");


	// Recover the primitive variables

	double rho_L = Lstate(0);
	double u_L = Lstate(1)/Lstate(0);
	double p_L = eos->p(Lstate);

	double rho_R = Rstate(0);
	double u_R = Rstate(1)/Rstate(0);
	double p_R = eos->p(Rstate);


	// Estimate pressure using PVRS scheme (linearising around arithmetic average state)

	double a_L = eos->a(Lstate);
	double a_R = eos->a(Rstate);
	double p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R-u_L)*0.5*(rho_L+rho_R)*0.5*(a_L+a_R);
	double p_star = std::max(0.0,p_pvrs);


	// Pressure-based wave speed estimates

	double S_L = u_L - a_L*q_K(eos->get_gamma(), p_star, p_L);
	double S_R = u_R + a_R*q_K(eos->get_gamma(), p_star, p_R);
	double S_star = (p_R-p_L+rho_L*u_L*(S_L-u_L)-rho_R*u_R*(S_R-u_R))/(rho_L*(S_L-u_L)-rho_R*(S_R-u_R));


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


double HLLC_RS_idealgas :: q_K (double gamma, double p_star, double p_K)
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






void exact_RS_idealgas :: solve_rp (	

	blitz::Array<double,1> Lstate,
	blitz::Array<double,1> Rstate,
	blitz::Array<double,1> flux,
	std::shared_ptr<eos_base> eos
)
{
	/*
	 * 	An exact Riemann solver for ideal gas. Flux across x/t=0 characteristic returned.
	 */

	assert(eos->get_eos_type() == "ideal");
	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));

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

	flux = euler_flux(soln(0), soln(1), soln(2), E);
}
	





// ========================================= RIEMANN SOLVERS FOR MIXED RIEMANN PROBLEM =====================================/


void M_HLLC_RS :: solve_rp_forinterfaceboundary (	

	blitz::Array<double,1> Lstate,
	blitz::Array<double,1> Rstate,
	double& p_star,
	double& u_star,
	double& rho_star_L,
	double& rho_star_R,
	std::shared_ptr<eos_base> eosL,
	std::shared_ptr<eos_base> eosR 
)
{
	/*
	 *	The M-HLLC solver of Hu and Adams (2009). An extension of HLLC to arbitrary EOS.
	 *	Finding unsatisfactory results with the star state density estimate, I have instead
	 *	computed densities using exact expressions for density change across waves.
	 */

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
	

	// The following approach from the paper gives unsatisfactory results!
	/*
	rho_star_L = rho_L*(u_L - S_L)/(u_star - S_L);
	
	rho_star_R = rho_R*(u_R - S_R)/(u_star - S_R);
	*/


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



double M_HLLC_RS :: mu (double fL, double fR, double rhoL, double rhoR)
{
	/*
	 *	Averaging used for many quantities so that U-property is satisfied
	 */

	return (sqrt(rhoL)*fL + sqrt(rhoR)*fR)/(sqrt(rhoL)+sqrt(rhoR));
}








void exact_RS_multi_idealgas :: solve_rp_forinterfaceboundary (	

	blitz::Array<double,1> Lstate,
	blitz::Array<double,1> Rstate,
	double& p_star,
	double& u_star,
	double& rho_star_L,
	double& rho_star_R,
	std::shared_ptr<eos_base> eosL,
	std::shared_ptr<eos_base> eosR 
)
{
	/*
	 * 	Exact mixed Riemann solver, where both left and right fluids obey the
	 *	ideal gas equation of state.
	 */

	assert(eosL->get_eos_type() == "ideal");
	assert(eosR->get_eos_type() == "ideal");

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
}
