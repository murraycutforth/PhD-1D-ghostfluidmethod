# The ghost fluid method for two-material flow governed by the one dimensional compressible Euler equations

*Murray Cutforth, Scientific Computing Group, University of Cambridge*



## Getting Started

These instructions will get you a copy of the project compiled on your local machine.

### Prerequisites

* The Blitz++ library - available at https://sourceforge.net/projects/blitz/
* Gnuplot - run `sudo apt-get install gnuplot`

### Installing

Choose your directory and run:

  `git clone https://github.com/murraycutforth/exact_riemann_solver_idealgas.git`  
  `git clone https://github.com/murraycutforth/1D_Euler_GFM.git`  
  `cd 1D_Euler_GFM`  
  `make`  

## Running a simulation

Simulation options are specified in 'settings_file.txt' before running 1D_Euler_GFM.exe. The following
options are available:

* length [integer] *number of real cells in domain*
* numGC [integer] *number of ghost cells on each end of the domain*
* lsnumGC [integer] *number of boundary ghost cells for level set grid*
* lslength [integer] *number of grid points for level set grid*
* fluid1_gamma [double] *EOS parameter*
* fluid1_A [double] *EOS parameter*
* fluid2_gamma [double] *EOS parameter*
* fluid2_A [double] *EOS parameter}
* RS_pure [HLLC_idealgas | Exact_idealgas] *Riemann solver used for single-material Riemann problem in flow solver*
* FS [Godunov | MUSCL] *Flow solver for single-material Euler equations*
* GFM [Original | Isobaricfix | Real] *Ghost fluid method*
* IC [TTC1 | TTC2 | TTC3 | TTC4 | TTC5] *Test case. TTC1-5 are initial conditions for the five test cases from Toro.*
* eos 1 [ideal] *Equation of state*
* eos 2 [ideal] *Equation of state*
* BC_L [transmissive | reflective | nothing] *Left boundary condition*
* BC_R [transmissive | reflective | nothing] *Right boundary condition*
* CFL [double between 0 and 1] *CFL number for time step size*
* sim [onefluid | twofluid] *Switch between single fluid simulation and full multimaterial ghost fluid simulation*
* outputpath [string] *Simulation output is sent here*

## Visualising the results

## Citations

* Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media, 2013.
* Fedkiw, Ronald P., et al. "A non-oscillatory Eulerian approach to interfaces in multimaterial flows (the ghost fluid method)." Journal of computational physics 152.2 (1999): 457-492.
* Hu, X. Y., N. A. Adams, and Gianluca Iaccarino. "On the HLLC Riemann solver for interface interaction in compressible multi-fluid flow." Journal of Computational Physics 228.17 (2009): 6572-6589.

