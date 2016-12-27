/*
 *	DESCRIPTION:	Data storage objects used in the simulation. The dimensions and position of an
 *			array is specified by an arrayinfo object. This finite volume method stores the
 *			integral-averaged vector of conserved quantities in each cell (the density of
 *			mass, momentum, and energy). Each fluid_state_array object also stores a pointer
 *			to an object describing the equation of state of this fluid. The level set_array
 *			represents a finite difference grid, where the value of the level set field is
 *			stored at every grid point.
 *
 *	TODO:		(Medium priority) Implement periodic boundary conditions in both classes - DONE (needs testing)
 *			(High priority) Implement linear extrapolation of level set in transmissive and reflective BCs
 *
 */


#include "data_storage.hpp"
#include "eos.hpp"
#include "misc.hpp"
#include <fstream>
#include <limits>
#include <cmath>
#include <cassert>


#define all blitz::Range::all()



double arrayinfo :: cellcentre_coord (int i)
{
	assert(i>=0);
	assert(i<=length + 2*numGC);

	return x0 - numGC*dx + 0.5*dx + i*dx;
}


int arrayinfo :: cellindex (double x)
{
	assert(x >= x0);
	assert(x <= x0 + length*dx);

	return static_cast<int>((x - x0)/dx) + numGC;
}


bool operator==(const arrayinfo& rhs, const arrayinfo& lhs)
{
	return (rhs.length == lhs.length)
		&& (rhs.x0 == lhs.x0)
		&& (rhs.dx == lhs.dx)
		&& (rhs.numGC == lhs.numGC);
}





fluid_state_array :: fluid_state_array (arrayinfo array, std::shared_ptr<eos_base> eos)
:
	array 	(array),
	CV 	(array.length + array.numGC*2, 3),
	eos 	(eos)
{}


fluid_state_array :: fluid_state_array (const fluid_state_array& other)
:
	array 	(other.array),
	CV	(other.CV),
	eos 	(other.eos)
{}


fluid_state_array fluid_state_array :: copy ()
{
	/*
	 *	Return a copy of this object with separate storage in memory
	 *	(since the copy constructor has reference semantics)
	 */
	
	fluid_state_array newstatearr (array, eos);
	newstatearr.CV = CV;
	return newstatearr;
}


void fluid_state_array :: apply_BCs()
{
	/*
	 * Set ghost cells on each end of array
	 */

	for (int i=0; i<array.numGC; i++)
	{
		if (array.leftBC == "transmissive")
		{
			CV(i,all) = CV(2*array.numGC - 1 - i,all);
		}
		else if (array.leftBC == "reflective")
		{
			CV(i,all) = CV(2*array.numGC - 1 - i,all);
			CV(i,1) = -CV(i,1);
		}
		else if (array.leftBC == "periodic")
		{
			CV(i,all) = CV(i+array.length,all);
		}
		else if (array.leftBC != "nothing")
		{
			assert(!"Invalid BC type");
		}
	}
	
	for (int i=0; i<array.numGC; i++)
	{
		if (array.rightBC == "transmissive")
		{
			CV(array.length + 2*array.numGC - 1 - i,all) = CV(array.length + i,all);
		}
		else if (array.rightBC == "reflective")
		{
			CV(array.length + 2*array.numGC - 1 - i,all) = CV(array.length + i,all);
			CV(array.length + 2*array.numGC - 1 - i,1) = -CV(array.length + 2*array.numGC - 1 - i,1);
		}
		else if (array.rightBC == "periodic")
		{
			CV(array.length + 2*array.numGC - 1 - i,all) = CV(2*array.numGC - 1 - i,all);
		}
		else if (array.rightBC != "nothing")
		{
			assert(!"Invalid BC type");
		}
	}
}


void fluid_state_array :: output_to_file (std::string name)
{
	/*
	 *	Output fluid density, velocity, total energy, pressure, and specific internal energy
	 */

	std::ofstream outfile;
	outfile.open(name);

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = array.cellcentre_coord(i);
		outfile << x << " " << CV(i,0) << " " << CV(i,1)/CV(i,0) << " " 
			<< CV(i,2) << " " << eos->p(CV(i,all))
			<< " " << specific_ie_cv(CV(i,all)) << std::endl;
	}

	outfile.close();
}


blitz::Array<double,1> fluid_state_array :: total_conserved_quantities ()
{
	/*
	 *	Return the integral of the vector of conserved quatities over all cells
	 */

	blitz::Array<double,1> sum (3);
	sum = 0.0;

	for (int i=array.numGC; i<array.length + array.numGC; i++)
	{
		sum(all) += array.dx*CV(i,all);
	}

	return sum;
}





levelset_array :: levelset_array (arrayinfo array)
:
	array	(array),
	phi	(array.length + 2*array.numGC)
{}


double levelset_array :: linear_interpolation (double x)
{
	/*
	 *	Interpolate to estimate level set value between grid points
	 */

	assert(array.numGC >= 1);
	assert(x >= array.cellcentre_coord(0));
	assert(x <= array.cellcentre_coord(array.length+2*array.numGC-1));


	// Find cell indices on L and R

	int i_L = static_cast<int>(floor((x - (array.x0 - 0.5*array.dx))/array.dx)) + array.numGC - 1;
	int i_R = i_L + 1;
	assert(i_L >= 0);
	assert(i_R <= array.length + 2*array.numGC-1);


	// Linear interpolation between values

	double t = x - array.cellcentre_coord(i_L);
	t = (t < 0.0) ? 0.0 : t;
	t = (t > array.dx) ? array.dx : t;

	return phi(i_L) + (t/array.dx)*(phi(i_R) - phi(i_L));
}


void levelset_array :: apply_BCs ()
{
	/*
	 * Set value of ghost cells on either end of array
	 */

	for (int i=0; i<array.numGC; i++)
	{
		if (array.leftBC == "transmissive")
		{
			double diff = phi(array.numGC+1) - phi(array.numGC);
			phi(i) = phi(array.numGC) - (array.numGC - i)*diff;
		}
		else if (array.leftBC == "reflective")
		{
			phi(i) = phi(2*array.numGC - 1 - i);
		}
		else if (array.leftBC == "nothing")
		{
			assert(!"Invalid BC type");
		}
	}
	
	for (int i=0; i<array.numGC; i++)
	{
		if (array.rightBC == "transmissive")
		{
			double diff = phi(array.length + array.numGC - 1) - phi(array.length + array.numGC - 2);
			phi(array.length + 2*array.numGC - 1 - i) = phi(array.length + array.numGC - 1) + (i - (array.length + array.numGC -1))*diff;
		}
		else if (array.rightBC == "reflective")
		{
			phi(array.length + 2*array.numGC - 1 - i) = phi(array.length + i);
		}
		else if (array.rightBC != "nothing")
		{
			assert(!"Invalid BC type");
		}
	}
}



void levelset_array :: output_to_file (std::string name)
{
	/*
	 *	Output level set grid to file
	 */

	std::ofstream outfile;
	outfile.open(name);

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = array.cellcentre_coord(i);
		outfile << x << " " << phi(i) << std::endl;
	}

	outfile.close();
}
