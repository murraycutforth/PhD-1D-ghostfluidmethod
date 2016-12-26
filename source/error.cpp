/*
 *	DESCRIPTION:	Functions used to calculate the error at the end of a simulation
 *
 */


#include "error.hpp"
#include "misc.hpp"
#include "exact_RS_idealgas.hpp"
#include <cassert>
#include <fstream>
#include <string>



blitz::Array<double,2> get_cellwise_error (
	
	fluid_state_array& fluid1,
	settingsfile& SF
)
{
	/*
	 *	If possible, this function computes the exact solution to the test
	 *	problem, and returns an array with the L1 error norm of the primitive
	 *	variables in each cell.
	 */

	blitz::Array<double,2> cellwise_error (fluid1.array.length,3);
	blitz::Array<double,1> leftprimitives (3);
	blitz::Array<double,1> rightprimitives (3);
	blitz::Array<double,1> soln (3);
	double discontinuitylocation;
	
	if (SF.IC == "TTC1")
	{
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1.0;
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.1;
		discontinuitylocation = 0.5;
	}
	else if (SF.IC == "TTC2")
	{
		leftprimitives(0) = 1.0;
		leftprimitives(1) = -2.0;
		leftprimitives(2) = 0.4;
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 2.0;
		rightprimitives(2) = 0.4;
		discontinuitylocation = 0.5;
	}
	else if (SF.IC == "TTC3")
	{
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1000.0;
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.01;
		discontinuitylocation = 0.5;
	}
	else if (SF.IC == "TTC4")
	{
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 0.01;
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 100.0;
		discontinuitylocation = 0.5;
	}
	else if (SF.IC == "TTC5")
	{
		leftprimitives(0) = 5.99924;
		leftprimitives(1) = 19.5975;
		leftprimitives(2) = 460.894;
		rightprimitives(0) = 5.99242;
		rightprimitives(1) = -6.19633;
		rightprimitives(2) = 46.0950;
		discontinuitylocation = 0.5;
	}
	else
	{
		assert(!"Invalid IC in error function");
	}
	

	exact_rs_idealgas RS (fluid1.eos->get_gamma(), fluid1.eos->get_gamma());
	RS.solve_RP(leftprimitives,rightprimitives);
	
	for (int i=0; i<fluid1.array.length; i++)
	{
		int materialcellind = i + fluid1.array.numGC;
		double x = fluid1.array.cellcentre_coord(materialcellind);
		double xot = (x - discontinuitylocation)/SF.T;
		soln = RS.sample_solution(xot);

		cellwise_error(i,0) = fabs(soln(0) - fluid1.CV(materialcellind,0));
		cellwise_error(i,1) = fabs(soln(1) - (fluid1.CV(materialcellind,1)/fluid1.CV(materialcellind,0)));
		cellwise_error(i,2) = fabs(soln(2) - fluid1.eos->p(fluid1.CV(materialcellind,blitz::Range::all())));
	}
	
	return cellwise_error;
}




void get_density_errornorms (

	blitz::Array<double,2> cellwise_error,
	double& L1error,
	double& Linferror
)
{
	/*
	 *	Using the array of errors of primitive variables in each cell, compute the L1 and L-infinity
	 *	error norms of the density.
	 */

	double sumerr = 0.0;
	double maxerr = 0.0;
	int N = cellwise_error.extent(blitz::firstDim);
	
	for (int i=0; i<N; i++)
	{
		sumerr += cellwise_error(i,0);	
		maxerr = std::max(maxerr, cellwise_error(i,0));
	}
	
	L1error = sumerr/N;
	Linferror = maxerr;
}



void output_errornorms_to_file (

	fluid_state_array& fluid1,
	settingsfile& SF
)
{
	/*
	 *	Store the L1 and Linf error in one file
	 */
	
	std::ofstream outfile;
	outfile.open(SF.basename + "finalerror.dat");

	blitz::Array<double,2> cellwise_error (get_cellwise_error(fluid1,SF));
	double L1err, Linferr;	
	get_density_errornorms(cellwise_error, L1err, Linferr);
	
	outfile << SF.length << " " << L1err << " " << Linferr << std::endl;
}
