#include "data_storage.hpp"
#include "eos.hpp"



#include <sstream>
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












void settingsfile :: read_settings_file ()
{
	std::ifstream infile("settings_file.txt");


	std::string BCRchoice;
		std::string RSchoice;
		std::string ICchoice;
		std::string GFMchoice;
		std::string simchoice;
		std::string eos1choice;
		std::string eos2choice;
		std::string FSchoice;



	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string inputname;

		iss >> inputname;

		if (inputname == "length") iss >> length;

		if (inputname == "x0") iss >> x0;

		if (inputname == "dx") iss >> dx;

		if (inputname == "numGC") iss >> numGC;

		if (inputname == "lsnumGC") iss >> lsnumGC;

		if (inputname == "lslength") iss >> lslength;

		if (inputname == "lsdx") iss >> lsdx;

		if (inputname == "fluid1_gamma") iss >> fluid1_gamma;

		if (inputname == "fluid1_B") iss >> fluid1_B;

		if (inputname == "fluid2_gamma") iss >> fluid2_gamma;

		if (inputname == "fluid2_B") iss >> fluid2_B;

		if (inputname == "T") iss >> T;

		if (inputname == "CFL") iss >> CFL;
		
		if (inputname == "RS")
		{
			iss >> RSchoice;
			if (RSchoice == "HLLC") RS = HLLC;
			else if (RSchoice == "M_HLLC") RS = M_HLLC;
			else if (RSchoice == "exact_idealgas") RS = exact_idealgas;
			else assert(!"Invalid RS");
		}
	
		if (inputname == "FS")
		{
			iss >> FSchoice;
			if (FSchoice == "Godunov") FS = Godunov;
			if (FSchoice == "MUSCL_FS") FS = MUSCL_FS;
		}

		if (inputname == "GFM")
		{
			iss >> GFMchoice;
			if (GFMchoice == "Original") GFM = Original;
			else if (GFMchoice == "Isobaricfix") GFM = Isobaricfix;
			else if (GFMchoice == "Real") GFM = Real;
			else assert(!"Invalid GFM");
		}

		if (inputname == "IC")
		{
			iss >> ICchoice;
			if (ICchoice == "TC1") IC = TC1;
			else if (ICchoice == "TC2") IC = TC2;
			else if (ICchoice == "HuST2") IC = HuST2;
			else if (ICchoice == "rGFMTC1") IC = rGFMTC1;
			else if (ICchoice == "rGFMTC3") IC = rGFMTC3;
			else assert(!"Invalid IC");
		}

		std::string lsICchoice;
		if (inputname == "lsIC")
		{
			iss >> lsICchoice;
			if (lsICchoice == "T1") lsIC = T1;
			else if (lsICchoice == "T2") lsIC = T2;
			else assert(!"Invalid lsIC");
		}

		if (inputname == "eos1")
		{
			iss >> eos1choice;
			if (eos1choice == "ideal") eos1 = ideal;
		}

		if (inputname == "eos2")
		{
			iss >> eos2choice;
			if (eos2choice == "ideal") eos2 = ideal;
		}

		std::string BCLchoice;
		if (inputname == "BC_L")
		{
			iss >> BCLchoice;
			if (BCLchoice == "reflective") BC_L = reflective;
			if (BCLchoice == "transmissive") BC_L = transmissive;
		}
		
		if (inputname == "BC_R")
		{
			iss >> BCRchoice;
			if (BCRchoice == "reflective") BC_R = reflective;
			if (BCRchoice == "transmissive") BC_R = transmissive;
		}

		if (inputname == "sim")
		{
			iss >> simchoice;
			if (simchoice == "serial_onefluid") sim = serial_onefluid;
			if (simchoice == "serial_twofluid") sim = serial_twofluid;
		}

		if (inputname == "outputpath") iss >> outputpath;

		

	}

	basename = outputpath + GFMchoice + "_" + FSchoice + "_" + ICchoice + "_" 
		+ eos1choice + "-" + eos2choice + "_" + std::to_string(length) + "_";

	if (simchoice == "serial_onefluid") basename = outputpath + "onefluid_" + FSchoice + "_" + ICchoice + "_" + eos1choice + "_" + std::to_string(length) + "_";
	
	infile.close();
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
	eos2	(eos2),
	initial_total_cv	(3)
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





void twofluid_array :: output_conservation_error (std::string name, levelset_array& ls, double t)
{
	std::ofstream outfile;

	if (t == 0.0)
	{
		outfile.open(name);
	}
	else
	{
		outfile.open(name, std::ios_base::app);
	}

	blitz::Array<double,1> total_cv = total_conserved_quantities(ls);

	outfile << t << " " << - initial_total_cv(0) + total_cv(0)
		<< " " << - initial_total_cv(1) + total_cv(1)
		<< " " << - initial_total_cv(2) + total_cv(2) << std::endl;

	outfile.close();
}





blitz::Array<double,1> twofluid_array :: total_conserved_quantities (levelset_array& ls)
{
	blitz::Array<double,1> sum (3);
	sum = 0.0;

	for (int i=array.numGC; i<array.length + array.numGC; i++)
	{
		if (ls(array.cellcentre_coord(i)) <= 0.0)
		{
			sum(all) += array.dx*fluid1(i,all);
		}
		else
		{
			sum(all) += array.dx*fluid2(i,all);
		}
	}

	return sum;
}



























onefluid_array :: onefluid_array (arrayinfo array, std::shared_ptr<eos_base> eos)
:
	fluid (array.length + array.numGC*2, 3),
	array (array),
	eos (eos),
	initial_total_cv	(3)
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
	sum = 0.0;

	for (int i=array.numGC; i<array.length + array.numGC; i++)
	{
		sum(all) += array.dx*fluid(i,all);
	}

	return sum;
}




void onefluid_array :: output_conservation_error (std::string name, double t)
{

	std::ofstream outfile;

	if (t == 0.0)
	{
		outfile.open(name);
	}
	else
	{
		outfile.open(name, std::ios_base::app);
	}

	blitz::Array<double,1> total_cv = total_conserved_quantities();

	outfile << t << " " << total_cv(0) - initial_total_cv(0)
		<< " " << total_cv(1) - initial_total_cv(1)
		<< " " << total_cv(2) - initial_total_cv(2) << std::endl;

	outfile.close();
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

	int i_L = static_cast<int>(floor((x - (array.x0 - 0.5*array.dx))/array.dx)) + array.numGC - 1;
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
