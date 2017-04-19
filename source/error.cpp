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
	else if (SF.IC == "GDA")
	{
		
		double u = 1.0;
		double p = 0.0001;
		double mu = 0.5;
		double A = 1000.0;
		double sigma = 0.1;

		for (int i=0; i<fluid1.array.length; i++)
		{
			int fluidcellind = i + fluid1.array.numGC;
			double x = fluid1.array.cellcentre_coord(fluidcellind);
			double exactrho = gaussian_function(A,mu,sigma,x);

			cellwise_error(i,0) = fabs(exactrho - fluid1.CV(fluidcellind,0));
			cellwise_error(i,1) = fabs(u - (fluid1.CV(fluidcellind,1)/fluid1.CV(fluidcellind,0)));
			cellwise_error(i,2) = fabs(p - fluid1.eos->p(fluid1.CV(fluidcellind,blitz::Range::all())));
		}

		return cellwise_error;		
	}
	else
	{
		assert(!"Invalid IC in error function");
	}
	

	exact_rs_idealgas RS (fluid1.eos->get_gamma(), fluid1.eos->get_gamma());
	RS.solve_RP(leftprimitives,rightprimitives);
	
	for (int i=0; i<fluid1.array.length; i++)
	{
		int fluidcellind = i + fluid1.array.numGC;
		double x = fluid1.array.cellcentre_coord(fluidcellind);
		double xot = (x - discontinuitylocation)/SF.T;
		soln = RS.sample_solution(leftprimitives, rightprimitives, xot);

		cellwise_error(i,0) = fabs(soln(0) - fluid1.CV(fluidcellind,0));
		cellwise_error(i,1) = fabs(soln(1) - (fluid1.CV(fluidcellind,1)/fluid1.CV(fluidcellind,0)));
		cellwise_error(i,2) = fabs(soln(2) - fluid1.eos->p(fluid1.CV(fluidcellind,blitz::Range::all())));
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



void get_velocity_errornorms (

	blitz::Array<double,2> cellwise_error,
	double& L1error,
	double& Linferror
)
{
	double sumerr = 0.0;
	double maxerr = 0.0;
	int N = cellwise_error.extent(blitz::firstDim);
	
	for (int i=0; i<N; i++)
	{
		sumerr += cellwise_error(i,1);	
		maxerr = std::max(maxerr, cellwise_error(i,1));
	}
	
	L1error = sumerr/N;
	Linferror = maxerr;
}



void get_pressure_errornorms (

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
		sumerr += cellwise_error(i,2);	
		maxerr = std::max(maxerr, cellwise_error(i,2));
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
	 *	Store the L1 and Linf error in density in one file
	 */
	
	

	blitz::Array<double,2> cellwise_error (get_cellwise_error(fluid1,SF));
	
	double L1errrho, Linferrrho;	
	get_density_errornorms(cellwise_error, L1errrho, Linferrrho);
	std::ofstream outfile;
	outfile.open(SF.basename + "densityerror.dat");
	outfile << SF.length << " " << L1errrho << std::endl;
	
	double L1erru, Linferru;	
	get_velocity_errornorms(cellwise_error, L1erru, Linferru);
	std::ofstream outfile2;
	outfile2.open(SF.basename + "velocityerror.dat");
	outfile2 << SF.length << " " << L1erru << std::endl;
	
	double L1errp, Linferrp;	
	get_pressure_errornorms(cellwise_error, L1errp, Linferrp);
	std::ofstream outfile3;
	outfile3.open(SF.basename + "pressureerror.dat");
	outfile3 << SF.length << " " << L1errp << std::endl;
}


void output_cellwise_error (

	fluid_state_array& fluid1,
	settingsfile& SF
)
{
	/*
	 *	Store the cellwise error in file
	 */

	std::ofstream outfile;
	outfile.open(SF.basename + "cellwiseerror.dat");

	blitz::Array<double,2> cellwise_error (get_cellwise_error(fluid1,SF));

	for (int i=0; i<SF.length; i++)
	{
		int fluidcellind = i + fluid1.array.numGC;
		double x = fluid1.array.cellcentre_coord(fluidcellind);
		outfile << x << " " << cellwise_error(i,0) << " " << cellwise_error(i,1) << " " << cellwise_error(i,2) << std::endl;
	}
}










void compute_total_U_onefluid (

	fluid_state_array& fluid1,
	blitz::Array<double,1> U0
)
{
	U0 = 0.0;

	for (int i=fluid1.array.numGC; i<fluid1.array.length + fluid1.array.numGC; i++)
	{
		U0 += fluid1.array.dx * fluid1.CV(i, blitz::Range::all());
	}
}


void update_total_U_onefluid (

	blitz::Array<double,1> FL,
	blitz::Array<double,1> FR,
	blitz::Array<double,1> U,
	double dt
)
{
	U += dt*(FL - FR);
}


void output_conservation_errors_to_file (
	
	blitz::Array<double,1> Ut,
	blitz::Array<double,1> U0,
	double t,
	settingsfile& SF
)
{
	std::ofstream outfile;

	if (t == 0.0) outfile.open(SF.basename + "conservationerror.dat");
	else outfile.open(SF.basename + "conservationerror.dat", std::fstream::app);

	blitz::Array<double,1> relative_error (3);
	relative_error(0) = (Ut(0) - U0(0))/U0(0);
	relative_error(1) = (Ut(1) - U0(1))/U0(1);
	relative_error(2) = (Ut(2) - U0(2))/U0(2);

	if (fabs(U0(1)) < 2e-16) relative_error(1) = 0.0;


	outfile << t << " " << relative_error(0) << " " << relative_error(1) << " " << relative_error(2) << std::endl;
}





void compute_total_U_twofluid (
	fluid_state_array& fluid1,
	fluid_state_array& fluid2,
	levelset_array& ls, 
	blitz::Array<double,1> U0
)
{
	U0 = 0.0;

	for (int i=fluid1.array.numGC; i<fluid1.array.length + fluid1.array.numGC; i++)
	{
		double fi = ls.linear_interpolation(fluid1.array.cellcentre_coord(i));

		if (fi <= 0.0)
		{
			U0 += fluid1.array.dx * fluid1.CV(i, blitz::Range::all());
		}
		else
		{
			U0 += fluid1.array.dx * fluid2.CV(i, blitz::Range::all());
		}
		
	}
}


void update_total_U_twofluid (

	blitz::Array<double,1> FL1,
	blitz::Array<double,1> FR1,
	blitz::Array<double,1> FL2,
	blitz::Array<double,1> FR2,
	levelset_array& ls,
	blitz::Array<double,1> U,
	double dt,
	fluid_state_array& fluid1
)
{
	int iL = fluid1.array.numGC;
	double phiL = ls.linear_interpolation(fluid1.array.cellcentre_coord(iL));

	int iR = fluid1.array.numGC + fluid1.array.length - 1;
	double phiR = ls.linear_interpolation(fluid1.array.cellcentre_coord(iR));

	if (phiL <= 0.0)
	{
		U += dt * FL1;
	}
	else
	{
		U += dt * FL2;
	}

	if (phiR <= 0.0)
	{
		U -= dt * FR1;
	}
	else
	{
		U -= dt * FR2;
	}
}