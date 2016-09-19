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
		rs->solve_rp(oldstate.fluid(i-1,all), oldstate.fluid(i,all), flux, edge_velocity(i-1), oldstate.eos, oldstate.eos);
		
		newstate.fluid(i-1, all) -= dtodx*flux;
		newstate.fluid(i, all) += dtodx*flux;
	}


}




