# The ghost fluid method for two-material flow governed by the one dimensional compressible Euler equations

Murray Cutforth, Scientific Computing Group, University of Cambridge



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

  length [integer above 0] {number of real cells in domain}
  numGC [integer greater than 0] {number of ghost cells on each end of the domain for boundary conditions}
  lsnumGC [integer greater than 0] {number of boundary ghost cells for level set grid}
  lslength [integer greater than 0] {number of grid points for level set grid}
  fluid1_gamma [double] {EOS parameter}
  fluid1_A [double] {EOS parameter}
  fluid2_gamma [double] {EOS parameter}
  fluid2_A [double] {EOS parameter}
  RS [HLLC | M_HLLC] {Riemann solver. HLLC for ideal gas using PVRS pressure based wave speeds, or M_HLLC for general EOS using  Roe-average state wave speed estimate}
  FS [Godunov | MUSCL_FS] {Flow solver. Godunov's first order method, or the second order MUSCL-Hancock method with minbee slope limiter}
  GFM [Original | Isobaricfix | Real] {Ghost fluid method. Original method, Original method with isobaric fix, real ghost fluid method [Wang, 2009]}
  IC [TC1 | TC2 | HuST2] {Test case. TC1 = Sod shock tube with 8x density ratio and 10x pressure ratio, x0=0.5. TC2 = two rarefaction test case, x0=0.5. HuST2 = Second shock tube test from Hu2009, 1xdensity ration, 2500x pressure ratio, gamma_R=1.667}
  eos 1 [ideal | tait]
  eos 2 [ideal | tait]
  BC_L [transmissive | reflective | nothing] {Left boundary condition}
  BC_R [transmissive | reflective | nothing] {Right boundary condition}
  CFL [double between 0 and 1] {CFL number for time step size}
  sim [onefluid | twofluid] {switches between different modes}
  outputpath [string] {simulation output is sent here}

## Visualising the results

## Citations

* Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media, 2013.
* Fedkiw, Ronald P., et al. "A non-oscillatory Eulerian approach to interfaces in multimaterial flows (the ghost fluid method)." Journal of computational physics 152.2 (1999): 457-492.
* Hu, X. Y., N. A. Adams, and Gianluca Iaccarino. "On the HLLC Riemann solver for interface interaction in compressible multi-fluid flow." Journal of Computational Physics 228.17 (2009): 6572-6589.

