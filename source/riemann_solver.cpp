#include "riemann_solver.hpp"
#include "flow_solver.hpp"
#include "eos.hpp"



#include <cassert>
#include <cmath>
#include <algorithm>






void HLLC_riemann_solver_idealgas ::  solve_rp (	blitz::Array<double,1> Lstate, 
							blitz::Array<double,1> Rstate, 
							blitz::Array<double,1> flux,
							double& S_star,
							std::shared_ptr<eos_base> eosL,
							std::shared_ptr<eos_base> eosR)
{
	assert(Lstate.extent(blitz::firstDim) == 3);
	assert(Rstate.extent(blitz::firstDim) == 3);
	assert(flux.extent(blitz::firstDim) == 3);
	assert(is_state_physical(Lstate));
	assert(is_state_physical(Rstate));
	assert(eosL->get_eos_type() == ideal);
	assert(eosR->get_eos_type() == ideal);
	assert(eosL->get_gamma() == eosL->get_gamma());
	assert(eosR->get_gamma() == eosR->get_gamma());


	// Recover the primitive variables

	double rho_L = Lstate(0);
	double u_L = Lstate(1)/Lstate(0);
	double p_L = eosL->p(Lstate);

	double rho_R = Rstate(0);
	double u_R = Rstate(1)/Rstate(0);
	double p_R = eosR->p(Rstate);


	// Estimate pressure using PVRS scheme (linearising around arithmetic average state)

	double a_L = eosL->a(Lstate);
	double a_R = eosR->a(Rstate);
	double p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R-u_L)*0.5*(rho_L+rho_R)*0.5*(a_L+a_R);
	double p_star = std::max(0.0,p_pvrs);


	// Pressure-based wave speed estimates

	double S_L = u_L - a_L*q_K(eosL->get_gamma(), p_star, p_L);
	double S_R = u_R + a_R*q_K(eosR->get_gamma(), p_star, p_R);
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
