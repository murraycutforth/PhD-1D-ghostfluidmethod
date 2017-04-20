/*
 *	DESCRIPTION:	Performs constant extrapolation away from the interface
 *			of the ghost states and the interface advection velocity.
 * 			The cells which are adjacent to the interface are held 
 * 			constant, and ghost cells are extrapolated outwards from
 * 			these values.
 *
 */


#include "extension_advection_eqn.hpp"
#include "misc.hpp"


#define all blitz::Range::all()


void extension_advection_eqn_1D (
	
	levelset_array& ls,
	fluid_state_array& ghoststatearr1,
	fluid_state_array& ghoststatearr2,
	blitz::Array<double,1> vfield,
	int frozenwidth
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

			if (!cell_local_to_interface(i, ghoststatearr1.array, ls, frozenwidth))
			{
				if (phi <= 0.0)
				{
					/*
					 *	Fluid 2 is ghost here
					 */

					if (normal > 0.0)
					{
						tempv = temp2.CV(i+1,all) - temp2.CV(i,all);
						delu = tempufield(i+1) - tempufield(i);
					}
					else
					{
						tempv = temp2.CV(i,all) - temp2.CV(i-1,all);
						delu = tempufield(i) - tempufield(i-1);
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
						tempv = temp1.CV(i,all) - temp1.CV(i-1,all);
						delu = tempufield(i) - tempufield(i-1);
					}
					else
					{
						tempv = temp1.CV(i+1,all) - temp1.CV(i,all);
						delu = tempufield(i+1) - tempufield(i);
					}

					ghoststatearr1.CV(i,all) -= dtodx*normal*tempv;
					vfield(i) -= dtodx*normal*delu;
				}
			}
		}
	}
}
