This file describes the options available in the settings file
length [integer above 0] {number of real cells in domain}
x0 [double] {coordinate of left edge of domain}
dx [double greater than 0] {width of each cell}
numGC [integer greater than 0] {number of ghost cells on each end of the domain for boundary conditions}
lsnumGC [integer greater than 0] {number of boundary ghost cells for level set grid}
lslength [integer greater than 0] {number of grid points for level set grid}
lsdx [double greater than 0] {length between level set grid points}
fluid1_gamma [double] {EOS parameter}
fluid1_A [double] {EOS parameter}
fluid2_gamma [double] {EOS parameter}
fluid2_A [double] {EOS parameter}
RS [HLLC | M_HLLC] {Riemann solver. HLLC for ideal gas using PVRS pressure based wave speeds, or M_HLLC for general EOS using Roe-average state wave speed estimate}
FS [Godunov | MUSCL_FS] {Flow solver. Godunov's first order method, or the second order MUSCL-Hancock method with minbee slope limiter}
GFM [Original | Isobaricfix] {Ghost fluid method. Original method, Original method with isobaric fix}
IC [TC1 | TC2 | HuST2] {Test case. TC1 = Sod shock tube with 8x density ratio and 10x pressure ratio, x0=0.5. TC2 = two rarefaction test case, x0=0.5. HuST2 = Second shock tube test from Hu2009, 1xdensity ration, 2500x pressure ratio, gamma_R=1.667}
lsIC [T1] {Initial conditions for level set. T1 places the zero contour at x=0.5, with fluid 1 to the right}
eos 1 [ideal]
eos 2 [ideal]
BC_L [transmissive | reflective | nothing]
BC_R [transmissive | reflective | nothing]
T [double greater than 0]
CFL [double between 0 and 1]
sim [serial_onefluid | serial_twofluid] {switches between different modes}
outputpath [string] {simulation output is sent here}
