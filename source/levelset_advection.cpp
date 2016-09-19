#include "levelset_advection.hpp"



#include <cassert>





void levelset_advection_firstorderup :: advection_operator (	levelset_array& ls,
								levelset_array& newls,
								std::shared_ptr<vfield_base> vfield,
								double dt)
{
	assert(dt > 0.0);
	assert(ls.array.numGC >= 1);
	assert(ls.array == newls.array);


	for (int i=ls.array.numGC; i<ls.array.length+ls.array.numGC; i++)
	{
		double cc = ls.array.cellcentre_coord(i);
		double u = vfield->get_u(cc);
		double phi_x;

		if (u > 0.0)
		{
			phi_x = (ls.phi(i) - ls.phi(i-1))/ls.array.dx;
		}
		else if (u < 0.0)
		{
			phi_x = (ls.phi(i+1) - ls.phi(i))/ls.array.dx;
		}
		else
		{
			phi_x = 0.0;
		}

		newls.phi(i) = ls.phi(i) - dt*u*phi_x;
	}
}

