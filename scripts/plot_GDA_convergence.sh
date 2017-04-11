# Create plot of error convergence for single fluid GDA test

IC=GDA
eos=ideal

gnuplot <<- EOF
	set terminal epslatex color size 9cm,9cm
	set palette rgb 33,13,10
	set size square
	set output "onefluid_${IC}_${eos}_errorconvergence.tex"
	set border lw 5
	set logscale x
	set logscale y
	set xrange [7:1500]
	set tics nomirror
	set xlabel 'Number of grid cells'
	set ylabel 'Density L1 error'
	plot 20000*x**(-1) w l ls 1 lw 4 title '\footnotesize First order convergence', 300000*x**(-2) w l ls 1 lc 2 lw 4 title '\footnotesize Second order convergence', "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 3 title '\footnotesize Godunov', "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 3 title '\footnotesize MUSCL-Hancock'

EOF
