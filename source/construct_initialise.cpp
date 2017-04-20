/*
 *	DESCRIPTION:	Functions used to construct and intiialise all objects needed to run a simulation
 *			in either single fluid or two fluid mode.
 *
 */

#include "construct_initialise.hpp"
#include "misc.hpp"
#include <cassert>
#include <string>

#define all blitz::Range::all()


fluid_state_array construct_initialise_onefluid ( 

	settingsfile& SF, 
	std::shared_ptr<eos_base>& eos, 
	std::shared_ptr<singlefluid_RS_base>& RS, 
	std::shared_ptr<flow_solver_base>& FS
)
{
	/*
	 *	Setup all simulation parameters using settings file, and return a
	 *	fully initialised material state array ready to begin time iterations
	 *	in a single fluid simulation.
	 */
	
	if (SF.eos1 == "ideal") eos = std::make_shared<eos_idealgas>(SF.fluid1_gamma);
	else assert(!"Invalid eos1");

	if (SF.RS_pure == "HLLC_idealgas") RS = std::make_shared<HLLC_RS_idealgas>();
	else if (SF.RS_pure == "Exact_idealgas") RS = std::make_shared<exact_RS_idealgas>();
	else assert(!"Invalid RS_pure");

	if (SF.FS == "Godunov") FS = std::make_shared<godunov>(RS);
	else if (SF.FS == "MUSCL") FS = std::make_shared<MUSCL>(RS);
	else assert(!"Invalid FS");

	arrayinfo array;
	array.length = SF.length;
	array.numGC = SF.numGC;
	array.leftBC = SF.BC_L;
	array.rightBC = SF.BC_R;


	// Use IC_type to set domain size
	
	if (SF.IC == "TTC1" || SF.IC == "TTC2" || SF.IC == "TTC3" || SF.IC == "TTC4" || SF.IC == "TTC5" || SF.IC == "GDA")
	{
		array.x0 = 0.0;
		array.dx = 1.0/array.length;
	}


	fluid_state_array statearr (array, eos);


	// Use IC_type to set initial values of fluid state

	if (SF.IC == "TTC1")
	{
		SF.T = 0.25;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.1;

		double discontinuitylocation = 0.5;

		set_piecewiseconstant_ICs ( leftprimitives, rightprimitives, discontinuitylocation, statearr);
	}
	else if (SF.IC == "TTC2")
	{
		SF.T = 0.15;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = -2.0;
		leftprimitives(2) = 0.4;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 2.0;
		rightprimitives(2) = 0.4;

		double discontinuitylocation = 0.5;

		set_piecewiseconstant_ICs ( leftprimitives, rightprimitives, discontinuitylocation, statearr);
	}
	else if (SF.IC == "TTC3")
	{
		SF.T = 0.012;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1000.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.01;

		double discontinuitylocation = 0.5;

		set_piecewiseconstant_ICs ( leftprimitives, rightprimitives, discontinuitylocation, statearr);
	}
	else if (SF.IC == "TTC4")
	{
		SF.T = 0.035;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 0.01;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 100.0;

		double discontinuitylocation = 0.5;

		set_piecewiseconstant_ICs ( leftprimitives, rightprimitives, discontinuitylocation, statearr);
	}
	else if (SF.IC == "TTC5")
	{
		SF.T = 0.035;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 5.99924;
		leftprimitives(1) = 19.5975;
		leftprimitives(2) = 460.894;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 5.99242;
		rightprimitives(1) = -6.19633;
		rightprimitives(2) = 46.0950;

		double discontinuitylocation = 0.5;

		set_piecewiseconstant_ICs ( leftprimitives, rightprimitives, discontinuitylocation, statearr);
	}
	else if (SF.IC == "GDA")
	{
		/*
		 *	GDA stands for Gaussian density advection.
		 */

		SF.T = 10.0;
		assert(SF.BC_R == "periodic");
		assert(SF.BC_L == "periodic");
		double u = 1.0;
		double p = 0.0001;
		double mu = 0.5;
		double A = 1000.0;
		double sigma = 0.1;

		for (int i=array.numGC; i<array.length+array.numGC; i++)
		{
			double x = array.cellcentre_coord(i);
			double rho = gaussian_function(A,mu,sigma,x);
			statearr.CV(i,0) = rho;
			statearr.CV(i,1) = rho*u;
			statearr.CV(i,2) = statearr.eos->E(rho,u,p);
		}
		
		statearr.apply_BCs();
	}
	else
	{
		assert(!"Invalid IC");
	}
	
	return statearr;
}



void set_piecewiseconstant_ICs (

	blitz::Array<double,1> leftprimitives,
	blitz::Array<double,1> rightprimitives,
	double discontinuitylocation,
	fluid_state_array& statearr
)
{
	/*
	 *	Given fluid primitive variables to the left and right, and the location
	 *	of the discontinuity between them, this function sets the conserved
	 *	variables in the statearr object.
	 */
	
	blitz::Array<double,1> leftstate (3);
	leftstate(0) = leftprimitives(0);
	leftstate(1) = leftprimitives(0)*leftprimitives(1);
	leftstate(2) = statearr.eos->E(leftprimitives);

	blitz::Array<double,1> rightstate (3);
	rightstate(0) = rightprimitives(0);
	rightstate(1) = rightprimitives(0)*rightprimitives(1);
	rightstate(2) = statearr.eos->E(rightprimitives);

	for (int i=statearr.array.numGC; i<statearr.array.length + statearr.array.numGC + 1; i++)
	{
		double x = statearr.array.cellcentre_coord(i);

		if (x < discontinuitylocation)
		{
			statearr.CV(i,all) = leftstate;
		}
		else
		{
			statearr.CV(i,all) = rightstate;
		}
	}

	statearr.apply_BCs();
}





void construct_initialise_twofluid (
	
	settingsfile& SF, 
	std::shared_ptr<eos_base>& eos1,
	std::shared_ptr<eos_base>& eos2, 
	std::shared_ptr<singlefluid_RS_base>& RS_pure,
	std::shared_ptr<multimat_RS_base>& RS_mixed, 
	std::shared_ptr<flow_solver_base>& FS,
	std::shared_ptr<GFM_base>& GFM,
	fluid_state_array& statearr1,
	fluid_state_array& statearr2,
	levelset_array& ls
)
{
	/*
	 *	Set up for a full ghost fluid method simulation
	 */
	
	if (SF.eos1 == "ideal") eos1 = std::make_shared<eos_idealgas>(SF.fluid1_gamma);
	else assert(!"Invalid eos1");
	if (SF.eos2 == "ideal") eos2 = std::make_shared<eos_idealgas>(SF.fluid2_gamma);
	else assert(!"Invalid eos2");
	
	if (SF.RS_pure == "HLLC_idealgas") RS_pure = std::make_shared<HLLC_RS_idealgas>();
	else if (SF.RS_pure == "Exact_idealgas") RS_pure = std::make_shared<exact_RS_idealgas>();
	else assert(!"Invalid RS_pure");
	if (SF.RS_mixed == "M_HLLC") RS_mixed = std::make_shared<M_HLLC_RS>();
	else if (SF.RS_mixed == "Exact_idealgas") RS_mixed = std::make_shared<exact_RS_multi_idealgas>();
	else assert(!"Invalid RS_mixed");

	if (SF.FS == "Godunov") FS = std::make_shared<godunov>(RS_pure);
	else if (SF.FS == "MUSCL") FS = std::make_shared<MUSCL>(RS_pure);
	else assert(!"Invalid FS");

	arrayinfo array;
	array.length = SF.length;
	array.numGC = SF.numGC;
	array.leftBC = SF.BC_L;
	array.rightBC = SF.BC_R;
	
	arrayinfo lsarray;
	lsarray.length = SF.lslength;
	lsarray.numGC = SF.lsnumGC;
	lsarray.leftBC = SF.BC_L;
	lsarray.rightBC = SF.BC_R;

	if (SF.IC == "TTC1" || SF.IC == "TTC2" || SF.IC == "TTC5")
	{
		array.x0 = 0.0;
		lsarray.x0 = 0.0;
		array.dx = 1.0/array.length;
		lsarray.dx = 1.0/lsarray.length;
	}

	statearr1.array = array;
	statearr1.eos = eos1;
	statearr1.CV.resize(array.length+2*array.numGC,3);
	statearr2.array = array;
	statearr2.eos = eos2;
	statearr2.CV.resize(array.length+2*array.numGC,3);	
	ls.array = lsarray;
	ls.phi.resize(lsarray.length+2*lsarray.numGC);

	if (SF.GFM == "OriginalGFM") GFM = std::make_shared<Original_GFM>(ls.array);
	else if (SF.GFM == "R_GFM") GFM = std::make_shared<R_GFM>(ls.array);
	else if (SF.GFM == "M_GFM") GFM = std::make_shared<M_GFM>(ls.array);
	else if (SF.GFM == "P_GFM") GFM = std::make_shared<P_GFM>(ls.array);
	else assert(!"Invalid GFM");

	if (SF.IC == "TTC1")
	{
		SF.T = 0.0007;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 100000.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 10000;

		double discontinuitylocation = 0.5;
		int parity = -1;

		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr1);
		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr2);
		set_singlediscontinuity_ls_IC (discontinuitylocation, parity, ls);
	}
	else if (SF.IC == "TTC2")
	{
		SF.T = 0.15;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = -2.0;
		leftprimitives(2) = 0.4;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 2.0;
		rightprimitives(2) = 0.4;

		double discontinuitylocation = 0.5;
		int parity = -1;

		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr1);
		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr2);
		set_singlediscontinuity_ls_IC (discontinuitylocation, parity, ls);
	}
	else if (SF.IC == "TTC5")
	{
		SF.T = 0.035;

		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 5.99924;
		leftprimitives(1) = 19.5975;
		leftprimitives(2) = 460.894;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 5.99242;
		rightprimitives(1) = -6.19633;
		rightprimitives(2) = 46.0950;
		
		double discontinuitylocation = 0.5;
		int parity = -1;

		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr1);
		set_piecewiseconstant_ICs (leftprimitives, rightprimitives, discontinuitylocation, statearr2);
		set_singlediscontinuity_ls_IC (discontinuitylocation, parity, ls);
	}
	else
	{
		assert(!"Invalid IC in SF");
	}
}
		



void set_singlediscontinuity_ls_IC (

	double discontinuitylocation, 
	int parity, 
	levelset_array& ls
)
{
	/*
	 *	Set up the level set with a single interface at the specified location.
	 *	The sign of the left hand region is given by parity.
	 */

	for (int i=ls.array.numGC; i<ls.array.numGC+ls.array.length; i++)
	{
		double x = ls.array.cellcentre_coord(i);
		double f = parity*(discontinuitylocation - x);
		ls.phi(i) = f;
	}

	ls.apply_BCs();
}		
