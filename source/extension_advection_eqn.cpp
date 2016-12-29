/*
 *	DESCRIPTION:	Performs constant extrapolation away from the interface
 *			of the ghost states and the interface advection velocity.
 *
 */


#include "extension_advection_eqn.hpp"
#include "misc.hpp"


#define all blitz::Range::all()


void extension_advection_eqn_1D (
	
	levelset_array& ls,
	fluid_state_array& ghoststatearr1,
	fluid_state_array& ghoststatearr2,
	blitz::Array<double,1> vfield
)
{
	int N = 10;
	double dtodx = 1.0;
	static fluid_state_array temp1 (ghoststatearr1.copy());
	static fluid_state_array temp2 (ghoststatearr2.copy());
	static blitz::Array<double,1> tempufield (vfield.extent(blitz::firstDim));
	blitz::Array<double,1> tempv (3);
	double delu;

	for (int n=0; n<N; n++)
	{
		temp1.CV = ghoststatearr1.CV;
		temp2.CV = ghoststatearr2.CV;
		tempufield = vfield;

		for (int i=ghoststatearr1.array.numGC; i<ghoststatearr1.array.length+ghoststatearr1.array.numGC; i++)
		{
			double phi = ls.linear_interpolation(ghoststatearr1.array.cellcentre_coord(i));
			double normal = ls.normal(ghoststatearr1.array.cellcentre_coord(i));

			if (!cell_local_to_interface(i, ghoststatearr1.array, ls))
			{
				if (phi <= 0.0)
				{
					/*
					 *	Fluid 2 is ghost here
					 */

					if (normal > 0.0)
					{
						tempv = ghoststatearr2.CV(i+1,all) - ghoststatearr2.CV(i,all);
						delu = vfield(i+1) - vfield(i);
					}
					else
					{
						tempv = ghoststatearr2.CV(i,all) - ghoststatearr2.CV(i-1,all);
						delu = vfield(i) - vfield(i-1);
					}

					ghoststatearr2.CV(i,all) += dtodx*normal*tempv;
					vfield(i) += dtodx*normal*delu;
				}
				else
				{
					/*
					 *	Fluid 1 is ghost here
					 */

					if (normal > 0.0)
					{
						tempv = ghoststatearr1.CV(i,all) - ghoststatearr1.CV(i-1,all);
						delu = vfield(i) - vfield(i-1);
					}
					else
					{
						tempv = ghoststatearr1.CV(i+1,all) - ghoststatearr1.CV(i,all);
						delu = vfield(i+1) - vfield(i);
					}

					ghoststatearr1.CV(i,all) -= dtodx*normal*tempv;
					vfield(i) -= dtodx*normal*delu;
				}
			}
		}
	}
}