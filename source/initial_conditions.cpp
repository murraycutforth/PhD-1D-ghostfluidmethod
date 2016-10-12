#include "initial_conditions.hpp"
#include "eos.hpp"



#include <cassert>
#include <cmath>



#define all blitz::Range::all()




void initialise_onefluid (onefluid_array& state, IC_type IC)
{
	if (IC == TC1)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.1;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = state.eos->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = state.eos->E(rightprimitives);

		double offsetstart = state.array.x0 - state.array.numGC*state.array.dx;

		for (int i=0; i<state.array.length + 2*state.array.numGC; i++)
		{
			double x = offsetstart + i*state.array.dx + 0.5*state.array.dx;

			if (x < 0.5)
			{
				state.fluid(i,all) = leftstate;
			}
			else
			{
				state.fluid(i,all) = rightstate;
			}
		}
	}
	else if (IC == TC2)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = -2.0;
		leftprimitives(2) = 0.4;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 2.0;
		rightprimitives(2) = 0.4;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = state.eos->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = state.eos->E(rightprimitives);

		double offsetstart = state.array.x0 - state.array.numGC*state.array.dx;

		for (int i=0; i<state.array.length + 2*state.array.numGC; i++)
		{
			double x = offsetstart + i*state.array.dx + 0.5*state.array.dx;

			if (x < 0.5)
			{
				state.fluid(i,all) = leftstate;
			}
			else
			{
				state.fluid(i,all) = rightstate;
			}
		}
	}
	else if (IC == HuST2)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 500;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.2;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = state.eos->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = state.eos->E(rightprimitives);

		double offsetstart = state.array.x0 - state.array.numGC*state.array.dx;

		for (int i=0; i<state.array.length + 2*state.array.numGC; i++)
		{
			double x = offsetstart + i*state.array.dx + 0.5*state.array.dx;

			if (x < 0.5)
			{
				state.fluid(i,all) = leftstate;
			}
			else
			{
				state.fluid(i,all) = rightstate;
			}
		}
	}
	else
	{
		assert(!"Invalid single fluid test case");
	}

	state.initial_total_cv = 0.0;

	for (int i=state.array.numGC; i<state.array.numGC+state.array.length; i++)
	{
		state.initial_total_cv(all) += state.array.dx*state.fluid(i,all);
	}
}









void initialise_twofluid (twofluid_array& states, IC_type IC)
{
	if (IC == TC1)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.1;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = states.eos2->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = states.eos1->E(rightprimitives);

		blitz::Array<double,1> zerostate (3);
		zerostate = 0.0;

		for (int i=0; i<states.array.length + 2*states.array.numGC; i++)
		{
			double x = states.array.cellcentre_coord(i);

			if (x < 0.5)
			{
				states.fluid2(i,all) = leftstate;
				states.fluid1(i,all) = zerostate;
			}
			else
			{
				states.fluid1(i,all) = rightstate;
				states.fluid2(i,all) = zerostate;
			}
		}
	}
	else if (IC == TC2)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = -2.0;
		leftprimitives(2) = 0.4;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 2.0;
		rightprimitives(2) = 0.4;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = states.eos2->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = states.eos1->E(rightprimitives);

		blitz::Array<double,1> zerostate (3);
		zerostate = 0.0;

		double offsetstart = states.array.x0 - states.array.numGC*states.array.dx;

		for (int i=0; i<states.array.length + 2*states.array.numGC; i++)
		{
			double x = offsetstart + i*states.array.dx + 0.5*states.array.dx;

			if (x < 0.5)
			{
				states.fluid2(i,all) = leftstate;
				states.fluid1(i,all) = zerostate;
			}
			else
			{
				states.fluid1(i,all) = rightstate;
				states.fluid2(i,all) = zerostate;
			}
		}
	}
	else if (IC == HuST2)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 500.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 0.2;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = states.eos2->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = states.eos1->E(rightprimitives);

		blitz::Array<double,1> zerostate (3);
		zerostate = 0.0;

		for (int i=0; i<states.array.length + 2*states.array.numGC; i++)
		{
			double x = states.array.cellcentre_coord(i);

			if (x < 0.5)
			{
				states.fluid2(i,all) = leftstate;
				states.fluid1(i,all) = zerostate;
			}
			else
			{
				states.fluid1(i,all) = rightstate;
				states.fluid2(i,all) = zerostate;
			}
		}
	}
	else if (IC == rGFMTC1)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 1.0;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 100000.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 0.125;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 10000.0;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = states.eos2->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = states.eos1->E(rightprimitives);

		blitz::Array<double,1> zerostate (3);
		zerostate = 0.0;

		for (int i=0; i<states.array.length + 2*states.array.numGC; i++)
		{
			double x = states.array.cellcentre_coord(i);

			if (x < 0.5)
			{
				states.fluid2(i,all) = leftstate;
				states.fluid1(i,all) = zerostate;
			}
			else
			{
				states.fluid1(i,all) = rightstate;
				states.fluid2(i,all) = zerostate;
			}
		}
	}
	else if (IC == rGFMTC3)
	{
		blitz::Array<double,1> leftprimitives (3);
		leftprimitives(0) = 0.82369077;
		leftprimitives(1) = 0.0;
		leftprimitives(2) = 1.0;

		blitz::Array<double,1> rightprimitives (3);
		rightprimitives(0) = 1.0;
		rightprimitives(1) = 0.0;
		rightprimitives(2) = 1.0;

		blitz::Array<double,1> leftstate (3);
		leftstate(0) = leftprimitives(0);
		leftstate(1) = leftprimitives(0)*leftprimitives(1);
		leftstate(2) = states.eos2->E(leftprimitives);

		blitz::Array<double,1> rightstate (3);
		rightstate(0) = rightprimitives(0);
		rightstate(1) = rightprimitives(0)*rightprimitives(1);
		rightstate(2) = states.eos1->E(rightprimitives);

		blitz::Array<double,1> zerostate (3);
		zerostate = 0.0;

		for (int i=0; i<states.array.length + 2*states.array.numGC; i++)
		{
			double x = states.array.cellcentre_coord(i);

			if (x < 0.2)
			{
				states.fluid2(i,all) = leftstate;
				states.fluid1(i,all) = zerostate;
			}
			else
			{
				states.fluid1(i,all) = rightstate;
				states.fluid2(i,all) = zerostate;
			}
		}
	}
	else
	{
		assert(!"Invalid test case");
	}


	states.initial_total_cv = 0.0;
	blitz::Array<double,1> zerostate (3);
	zerostate = 0.0;

	for (int i=states.array.numGC; i<states.array.numGC+states.array.length; i++)
	{
		if (states.fluid1(i,0) == 0.0 && states.fluid1(i,1) == 0.0 && states.fluid1(i,2) == 0.0)
		{
			states.initial_total_cv(all) += states.array.dx*states.fluid2(i,all);
		}
		else
		{
			states.initial_total_cv(all) += states.array.dx*states.fluid1(i,all);
		}
	}
}













void initialise_levelset (levelset_array& ls, ls_IC_type IC)
{
	if (IC == T1)
	{

		double offsetstart = ls.array.x0 - ls.array.numGC*ls.array.dx;

		for (int i=0; i<ls.array.length + 2*ls.array.numGC; i++)
		{
			double x = offsetstart + i*ls.array.dx + 0.5*ls.array.dx;

			if (x < 0.5)
			{
				ls.phi(i) = fabs(x - 0.5);
			}
			else
			{
				ls.phi(i) = -fabs(x - 0.5);
			}
		}
	}
	if (IC == T2)
	{

		double offsetstart = ls.array.x0 - ls.array.numGC*ls.array.dx;

		for (int i=0; i<ls.array.length + 2*ls.array.numGC; i++)
		{
			double x = offsetstart + i*ls.array.dx + 0.5*ls.array.dx;

			if (x < 0.2)
			{
				ls.phi(i) = fabs(x - 0.5);
			}
			else
			{
				ls.phi(i) = -fabs(x - 0.5);
			}
		}
	}
}

