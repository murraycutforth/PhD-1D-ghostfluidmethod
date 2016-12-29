# Given the IC name, create plot of error convergence

IC=$1
eos=ideal

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,35' size 2500,2500
	set output "./../output/onefluid_${IC}_${eos}_errorconvergence.png"
	set border lw 5
	set title "L1 error norm convergence of single fluid methods on test case: ${IC}" font 'Palatino,50'
	set logscale x
	set logscale y
	set tics nomirror
	set xlabel 'Number of grid cells'
	set xrange [7:1500]
	plot x**(-1) w l lw 3 title "First order convergence", "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collatedfinalerrors.dat" w p pt 7 ps 5 title "Godunov's method with exact riemann solver", "./../output/onefluid_Godunov_HLLC_idealgas_${IC}_${eos}_collatedfinalerrors.dat" w p pt 7 ps 5 title "Godunov's method with HLLC riemann solver", "./../output/onefluid_MUSCL_HLLC_idealgas_${IC}_${eos}_collatedfinalerrors.dat" w p pt 7 ps 5 title "MUSCL-Hancock method with HLLC riemann solver"

EOF
