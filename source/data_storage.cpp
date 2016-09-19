#include "data_storage.hpp"
#include "eos.hpp"



#include <fstream>
#include <limits>
#include <cmath>




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















twofluid_array :: twofluid_array (arrayinfo array, std::shared_ptr<eos_base> eos1, std::shared_ptr<eos_base> eos2)
:
	fluid1 (array.length + 2*array.numGC,3),
	fluid2 (array.length + 2*array.numGC,3),
	array	(array),
	eos1	(eos1),
	eos2	(eos2)
{}

	

void twofluid_array :: apply_BCs ()
{
	for (int i=0; i<array.numGC; i++)
	{
		if (array.leftBC == transmissive)
		{
			fluid1(i,all) = fluid1(2*array.numGC - 1 - i,all);
			fluid2(i,all) = fluid2(2*array.numGC - 1 - i,all);
		}
		else if (array.leftBC == reflective)
		{
			fluid1(i,all) = fluid1(2*array.numGC - 1 - i,all);
			fluid1(i,1) = -fluid1(i,1);
			fluid2(i,all) = fluid2(2*array.numGC - 1 - i,all);
			fluid2(i,1) = -fluid2(i,1);
		}
		else if (array.leftBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
	
	for (int i=0; i<array.numGC; i++)
	{
		if (array.rightBC == transmissive)
		{
			fluid1(array.length + 2*array.numGC - 1 - i,all) = fluid1(array.length + i,all);
			fluid2(array.length + 2*array.numGC - 1 - i,all) = fluid2(array.length + i,all);
		}
		else if (array.rightBC == reflective)
		{
			fluid1(array.length + 2*array.numGC - 1 - i,all) = fluid1(array.length + i,all);
			fluid1(array.length + 2*array.numGC - 1 - i,1) = -fluid1(array.length + 2*array.numGC - 1 - i,1);
			fluid2(array.length + 2*array.numGC - 1 - i,all) = fluid2(array.length + i,all);
			fluid2(array.length + 2*array.numGC - 1 - i,1) = -fluid2(array.length + 2*array.numGC - 1 - i,1);
		}
		else if (array.rightBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
}



void twofluid_array :: output_to_file (std::string name)
{
	std::ofstream outfile;
	outfile.open(name);

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = array.cellcentre_coord(i);
		outfile << x << " " << fluid1(i,0) << " " << fluid1(i,1)/fluid1(i,0) << " " 
			<< fluid1(i,2) << " " << eos1->p(fluid1(i,all))
			<< " " << specific_ie_cv(fluid1(i,all));
		outfile << " " << fluid2(i,0) << " " << fluid2(i,1)/fluid2(i,0) << " " 
			<< fluid2(i,2) << " " << eos2->p(fluid2(i,all))
			<< " " << specific_ie_cv(fluid2(i,all)) << std::endl;
	}

	outfile.close();
}





void twofluid_array :: output_realfluid_to_file (std::string name, levelset_array& ls)
{
	assert(ls.array.cellcentre_coord(0) < array.cellcentre_coord(0));
	assert(ls.array.cellcentre_coord(ls.array.length+2*ls.array.numGC-1)
		> array.cellcentre_coord(array.length+2*array.numGC-1));

	std::ofstream outfile;
	outfile.open(name);

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = array.cellcentre_coord(i);

		if (ls(x) <= 0.0)
		{
			outfile << x << " " << fluid1(i,0) << " " << fluid1(i,1)/fluid1(i,0) << " " 
				<< fluid1(i,2) << " " << eos1->p(fluid1(i,all))
				<< " " << specific_ie_cv(fluid1(i,all)) << std::endl;
		}
		else
		{
			outfile << x << " " << fluid2(i,0) << " " << fluid2(i,1)/fluid2(i,0) << " " 
				<< fluid2(i,2) << " " << eos2->p(fluid2(i,all))
				<< " " << specific_ie_cv(fluid2(i,all)) << std::endl;
		}
	}

	outfile.close();
}


















onefluid_array :: onefluid_array (arrayinfo array, std::shared_ptr<eos_base> eos)
:
	fluid (array.length + array.numGC*2, 3),
	array (array),
	eos (eos)
{}




onefluid_array :: onefluid_array (blitz::Array<double,2>& fluid, arrayinfo array, std::shared_ptr<eos_base> eos)
:
	fluid (fluid),
	array (array),
	eos (eos)
{}




void onefluid_array :: apply_BCs()
{
	for (int i=0; i<array.numGC; i++)
	{
		if (array.leftBC == transmissive)
		{
			fluid(i,all) = fluid(2*array.numGC - 1 - i,all);
		}
		else if (array.leftBC == reflective)
		{
			fluid(i,all) = fluid(2*array.numGC - 1 - i,all);
			fluid(i,1) = -fluid(i,1);
		}
		else if (array.leftBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
	
	for (int i=0; i<array.numGC; i++)
	{
		if (array.rightBC == transmissive)
		{
			fluid(array.length + 2*array.numGC - 1 - i,all) = fluid(array.length + i,all);
		}
		else if (array.rightBC == reflective)
		{
			fluid(array.length + 2*array.numGC - 1 - i,all) = fluid(array.length + i,all);
			fluid(array.length + 2*array.numGC - 1 - i,1) = -fluid(array.length + 2*array.numGC - 1 - i,1);
		}
		else if (array.rightBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
}





void onefluid_array :: output_to_file (std::string name)
{
	std::ofstream outfile;
	outfile.open(name);

	double offset = array.x0 - array.numGC*array.dx;

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = offset + i*array.dx + 0.5*array.dx;
		outfile << x << " " << fluid(i,0) << " " << fluid(i,1)/fluid(i,0) << " " 
			<< fluid(i,2) << " " << eos->p(fluid(i,all))
			<< " " << specific_ie_cv(fluid(i,all)) << std::endl;
	}

	outfile.close();
}




blitz::Array<double,1> onefluid_array :: total_conserved_quantities ()
{
	blitz::Array<double,1> sum (3);

	for (int i=array.numGC; i<array.length + array.numGC; i++)
	{
		sum(all) += array.dx*fluid(i,all);
	}

	return sum;
}

























levelset_array :: levelset_array (arrayinfo array)
:
	phi	(array.length + 2*array.numGC),
	array	(array)
{}



double levelset_array :: operator() (double x)
{
	assert(x >= array.cellcentre_coord(0));
	assert(x <= array.cellcentre_coord(array.length+2*array.numGC-1));

	return linear_interpolation(x);
}



double levelset_array :: linear_interpolation (double x)
{
	assert(array.numGC >= 1);


	// Find cell indices on L and R

	int i_L = static_cast<int>((x - (array.x0 - 0.5*array.dx))/array.dx) + array.numGC - 1;
	int i_R = i_L + 1;
	assert(i_L >= 0);
	assert(i_R <= array.length + 2*array.numGC-1);


	// Linear interpolation between values

	double t = x - array.cellcentre_coord(i_L);
	if (fabs(t) < 1e-10) t = 0.0;
	if (fabs(t - array.dx) < 1e-10) t = array.dx;
	assert(t >= 0.0);
	assert(t <= array.dx);

	return phi(i_L) + (t/array.dx)*(phi(i_R) - phi(i_L));
}




void levelset_array :: apply_BCs ()
{
	for (int i=0; i<array.numGC; i++)
	{
		if (array.leftBC == transmissive)
		{
			phi(i) = phi(array.numGC);
		}
		else if (array.leftBC == reflective)
		{
			phi(i) = phi(2*array.numGC - 1 - i);
		}
		else if (array.leftBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
	
	for (int i=0; i<array.numGC; i++)
	{
		if (array.rightBC == transmissive)
		{
			phi(array.length + 2*array.numGC - 1 - i) = phi(array.length + array.numGC - 1);
		}
		else if (array.rightBC == reflective)
		{
			phi(array.length + 2*array.numGC - 1 - i) = phi(array.length + i);
		}
		else if (array.rightBC == nothing)
		{
		}
		else
		{
			assert(!"Invalid BC type");
		}
	}
}



void levelset_array :: output_to_file (std::string name)
{
	std::ofstream outfile;
	outfile.open(name);

	double offset = array.x0 - array.numGC*array.dx;

	for (int i=0; i<array.length + 2*array.numGC; i++)
	{
		double x = offset + i*array.dx + 0.5*array.dx;
		outfile << x << " " << phi(i) << std::endl;
	}

	outfile.close();
}
