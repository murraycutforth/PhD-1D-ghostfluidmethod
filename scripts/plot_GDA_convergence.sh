# Create plot of error convergence for single fluid GDA test

IC=GDA
eos=ideal

gnuplot <<- EOF
	set terminal epslatex color size 16cm,8cm
	set output "onefluid_${IC}_${eos}_errorconvergence.tex"
	set multiplot layout 1,2
	set border lw 3
	set size square
	set tics nomirror

	
	set xrange [0:1]
	set xlabel 'x'
	set yrange [0:2000]
	set ylabel 'Density'
	plot 1000*exp(-((x-0.5)*(x-0.5))/(2*0.01)) w l ls 2 lw 2 lc 0 title "Exact", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_10_final.dat" u 1:2 w l ls 1 lw 2 lc 1 title "N = 10", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_40_final.dat" u 1:2 w l ls 1 lw 2 lc 2 title "N = 40", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_160_final.dat" u 1:2 w l ls 1 lw 2 lc 3 title "N = 160", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_640_final.dat" u 1:2 w l ls 1 lw 2 lc 4 title "N = 640"


	set logscale x
	set logscale y
	set xrange [7:1500]
	set yrange [0.1:50000]
	set ytics ("\$10^{-4}\$"0.0001,"\$10^{-3}\$"0.001, "\$10^{-2}\$" 0.01,"\$10^{-1}\$"0.1, "\$10^{0}\$" 1, "\$10^{1}\$" 10, "\$10^{2}\$" 100, "\$10^{3}\$" 1000, "\$10^{4}\$" 10000, "\$10^{5}\$" 100000)
	set xlabel 'Number of grid cells'
	set ylabel 'Density L1 error'
	plot 20000*x**(-1) w l ls 1 lw 4 title '\footnotesize First order convergence', 300000*x**(-2) w l ls 1 lc 2 lw 4 title '\footnotesize Second order convergence', "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 3 title '\footnotesize Godunov', "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 3 title '\footnotesize MUSCL-Hancock'

EOF
