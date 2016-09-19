#include "vfield.hpp"




vfield_base :: vfield_base ()
:
	t	(0.0)
{}











vfield_test1 :: vfield_test1 ()
:
	vfield_base()
{}





double vfield_test1 :: get_u (double x)
{
	return 1.0;
}











vfield_oldstate :: vfield_oldstate (twofluid_array& oldstates, levelset_array& oldls)
:
	states	(oldstates),
	ls	(oldls)
{}





double vfield_oldstate :: get_u (double x)
{
	// Interpolate linearly between real velocities

	int cellindex = states.array.cellindex(x);
	
	int indexL, indexR;

	if (x > states.array.cellcentre_coord(cellindex))
	{
		indexL = cellindex;
		indexR = cellindex + 1;
	}
	else
	{
		indexL = cellindex - 1;
		indexR = cellindex;
	}

	double u_L, u_R;

	if (ls(states.array.cellcentre_coord(indexL)) <= 0.0)
	{
		// fluid 1 is real at indexL

		u_L = states.fluid1(indexL,1)/states.fluid1(indexL,0);
	}
	else
	{
		u_L = states.fluid2(indexL,1)/states.fluid2(indexL,0);
	}

	if (ls(states.array.cellcentre_coord(indexR)) <= 0.0)
	{
		u_R = states.fluid1(indexR,1)/states.fluid1(indexR,0);
	}
	else
	{
		u_R = states.fluid2(indexR,1)/states.fluid2(indexR,0);
	}

	double t = x - states.array.cellcentre_coord(indexL);
	if (fabs(t) < 1e-10) t = 0.0;
	if (fabs(t - states.array.dx) < 1e-10) t = states.array.dx;
	assert(t >= 0.0);
	assert(t <= states.array.dx);

	return u_L + (t/states.array.dx)*(u_R - u_L);
}











vfield_starstate :: vfield_starstate(	std::shared_ptr<flow_solver_base> FS,
				levelset_array& ls,
				arrayinfo statearray)
:
	FS	(FS),
	ls	(ls),
	statearray (statearray)
{}



double vfield_starstate :: get_u (double x)
{
	// Move along level set in both directions until a sign change is found

		// Find the closest edge in the edge_velocity array, return value
	
	int baseind = ls.array.cellindex(x);
	
	for (int i=1; i<ls.array.length+2*ls.array.numGC; i++)
	{
		int i_LL = baseind - i;
		int i_L = baseind - i + 1;
		int i_RR = baseind + i;
		int i_R = baseind + i - 1;

		if (0 <= i_LL && i_LL < ls.array.length + 2*ls.array.numGC
			&& 0 <= i_L && i_L < ls.array.length + 2*ls.array.numGC)
		{
			if (std::copysign(1.0,ls.phi(i_LL)) != std::copysign(1.0,ls.phi(i_L)))
			{
				double Linterface = ls.array.cellcentre_coord(i_LL);

				int Lindex = statearray.cellindex(Linterface);
				
				return FS->edge_velocity(Lindex);
			}
		}

		if (0 <= i_RR && i_RR < ls.array.length + 2*ls.array.numGC
			&& 0 <= i_R && i_R < ls.array.length + 2*ls.array.numGC)
		{
			if (std::copysign(1.0,ls.phi(i_RR)) != std::copysign(1.0,ls.phi(i_R)))
			{

				double Linterface = ls.array.cellcentre_coord(i_R);

				int Lindex = statearray.cellindex(Linterface);
				
				return FS->edge_velocity(Lindex);
			}
		}
	}

	assert(false);
}
