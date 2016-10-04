#include "flow_solver.hpp"
#include "eos.hpp"
#include "riemann_solver.hpp"




#define all blitz::Range::all()






blitz::Array<double,1> euler_flux (double rho, double u, double P, double E)
{
	blitz::Array<double,1> flux (3);

	flux(0) = rho*u;
	flux(1) = rho*u*u + P;
	flux(2) = u*(E+P);

	return flux;
}
	




blitz::Array<double,1> euler_flux (blitz::Array<double,1> cv, std::shared_ptr<eos_base> eos)
{

	blitz::Array<double,1> flux (3);

	double P = eos->p(cv);
	double u = cv(1)/cv(0); 

	flux(0) = cv(1);
	flux(1) = cv(1)*u + P;
	flux(2) = u*(cv(2)+P);

	return flux;
}











double limited_slope (double beta, double del_L, double del_R)
{
	double result;

	if (del_R > 0.0)
	{
		result = std::max(0.0, std::min(beta*del_L, del_R));
		result = std::max(result, std::min(del_L, beta*del_R));
	}
	else
	{
		result = std::min(0.0, std::max(beta*del_L, del_R));
		result = std::min(result, std::max(del_L, beta*del_R));
	}

	return result;
}




flow_solver_base :: flow_solver_base (std::shared_ptr<riemann_solver_base> rs)
:
	edge_velocity (),
	rs (rs)
{}




godunov :: godunov (std::shared_ptr<riemann_solver_base> rs)
:
	flow_solver_base (rs)
{}




void godunov :: single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt)
{
	assert(oldstate.fluid.extent(blitz::firstDim) == newstate.fluid.extent(blitz::firstDim));
	assert(oldstate.fluid.extent(blitz::secondDim) == newstate.fluid.extent(blitz::secondDim));
	assert(oldstate.fluid.extent(blitz::secondDim) == 3);
	assert(oldstate.array == newstate.array);
	assert(dt > 0.0);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 1);
	
	
	double dtodx = dt/oldstate.array.dx;
	blitz::Array<double,1> flux (3);


	// Ensure storage for edge velocities is correct size

	if (edge_velocity.extent(blitz::firstDim) != oldstate.array.length + 2*oldstate.array.numGC-1)
	{
		edge_velocity.resize(oldstate.array.length + 2*oldstate.array.numGC -1);
	}

	
	newstate.fluid = oldstate.fluid;






	// Update all real cells correctly

	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC + 1; i++)
	{
		rs->solve_rp(oldstate.fluid(i-1,all), oldstate.fluid(i,all), flux, edge_velocity(i-1), oldstate.eos);
		
		newstate.fluid(i-1, all) -= dtodx*flux;
		newstate.fluid(i, all) += dtodx*flux;

	}	
}

















MUSCL :: MUSCL (std::shared_ptr<riemann_solver_base> rs)
:
	flow_solver_base (rs)
{}



void MUSCL :: single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt)
{


	assert(oldstate.fluid.extent(blitz::firstDim) == newstate.fluid.extent(blitz::firstDim));
	assert(oldstate.fluid.extent(blitz::secondDim) == newstate.fluid.extent(blitz::secondDim));
	assert(oldstate.fluid.extent(blitz::secondDim) == 3);
	assert(oldstate.array == newstate.array);
	assert(dt > 0.0);
	assert(oldstate.array.length > 2*oldstate.array.numGC);
	assert(oldstate.array.numGC >= 2);
	
	
	double dtodx = dt/oldstate.array.dx;
	blitz::Array<double,1> flux (3);


	// Ensure storage for edge velocities is correct size

	if (edge_velocity.extent(blitz::firstDim) != oldstate.array.length + 2*oldstate.array.numGC-1)
	{
		edge_velocity.resize(oldstate.array.length + 2*oldstate.array.numGC -1);
	}

	
	blitz::Array<double,1> L_del_L (3);
	blitz::Array<double,1> L_del_R (3);
	blitz::Array<double,1> R_del_L (3);
	blitz::Array<double,1> R_del_R (3);
	blitz::Array<double,1> slope_L (3);
	blitz::Array<double,1> slope_R (3);
	blitz::Array<double,1> L_BEV_L (3);
	blitz::Array<double,1> L_BEV_R (3);
	blitz::Array<double,1> R_BEV_L (3);
	blitz::Array<double,1> R_BEV_R (3);
	blitz::Array<double,1> L_BEV_R_evolved (3);
	blitz::Array<double,1> R_BEV_L_evolved (3);

	

	newstate.fluid = oldstate.fluid;




	for (int i=oldstate.array.numGC; i<oldstate.array.length + oldstate.array.numGC + 1; i++)
	{
		// Left cell has index i-1, right cell has index i
		
		L_del_L = oldstate.fluid(i-1,all) - oldstate.fluid(i-2,all);
		L_del_R = oldstate.fluid(i,all) - oldstate.fluid(i-1,all);
		R_del_L = oldstate.fluid(i,all) - oldstate.fluid(i-1,all);
		R_del_R = oldstate.fluid(i+1,all) - oldstate.fluid(i,all);

		slope_L = MUSCL_slope (	0.5, 
					L_del_L,
					L_del_R);
		slope_R = MUSCL_slope (	0.5,
					R_del_L,
					R_del_R);


		// Limit slopes
		
		double beta = 1.0;
		for (int k=0; k<3; k++)
		{
			slope_L(k) = limited_slope(beta,L_del_L(k),L_del_R(k));
			slope_R(k) = limited_slope(beta,R_del_L(k),R_del_R(k));
		}
		

		// Construct boundary extrapolated values for both cells

		L_BEV_L = oldstate.fluid(i-1,all) - 0.5*slope_L;
		L_BEV_R = oldstate.fluid(i-1,all) + 0.5*slope_L;
		R_BEV_L = oldstate.fluid(i,all) - 0.5*slope_R;
		R_BEV_R = oldstate.fluid(i,all) + 0.5*slope_R;


		// Evolve states by half a time step

		L_BEV_R_evolved = L_BEV_R + 0.5*(dtodx)*(euler_flux(L_BEV_L, oldstate.eos) - euler_flux(L_BEV_R, oldstate.eos));

		R_BEV_L_evolved = R_BEV_L + 0.5*(dtodx)*(euler_flux(R_BEV_L, oldstate.eos) - euler_flux(R_BEV_R, oldstate.eos));


		// Update using conventional Riemann problem solution

		rs->solve_rp(L_BEV_R_evolved, R_BEV_L_evolved, flux, edge_velocity(i-1), oldstate.eos);
		
		newstate.fluid(i-1, all) -= dtodx*flux;
		newstate.fluid(i, all) += dtodx*flux;
	}
}



blitz::Array<double,1> MUSCL :: MUSCL_slope (double omega, blitz::Array<double,1> delU_L, blitz::Array<double,1> delU_R)
{
	blitz::Array<double,1> result (3);
	
	result = 0.5*(1.0+omega)*delU_L + 0.5*(1.0*omega)*delU_R;
	return result;
}













debugger :: debugger (std::shared_ptr<riemann_solver_base> rs)
:
	flow_solver_base(rs)
{}



void debugger :: single_fluid_update (onefluid_array& oldstate, onefluid_array& newstate, double dt)
{
}
