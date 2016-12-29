# Given the IC name, create plot of error convergence

IC=$1
GFM=$2

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,35' size 2500,2500
	set output "./../output/twofluid_${GFM}_${IC}_errorconvergence.png"
	set border lw 5
	set title "L1 error norm of GFM on test case: ${IC} using GFM: ${GFM}" font 'Palatino,50'
	set logscale x
	set logscale y
	set tics nomirror
	set xlabel 'Number of grid cells'
	set xrange [7:1500]
	plot x**(-1) w l lw 3 title "First order convergence", "./../output/${GFM}_Godunov_Exact_idealgas_Exact_idealgas_${IC}_ideal-ideal_collatedfinalerrors.dat" w p pt 7 ps 5 title "Godunov's method with exact riemann solver", "./../output/${GFM}_Godunov_HLLC_idealgas_M_HLLC_${IC}_ideal-ideal_collatedfinalerrors.dat" w p pt 7 ps 5 title "Godunov's method with HLLC riemann solver", "./../output/${GFM}_MUSCL_HLLC_idealgas_M_HLLC_${IC}_ideal-ideal_collatedfinalerrors.dat" w p pt 7 ps 5 title "MUSCL-Hancock method with HLLC riemann solver"

EOF
