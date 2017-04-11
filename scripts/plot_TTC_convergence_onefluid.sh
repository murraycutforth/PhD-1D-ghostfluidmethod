# Given the IC name, create plot of error convergence

IC=$1
eos=ideal

gnuplot <<- EOF
	set terminal epslatex color size 16cm,16cm
	set output "onefluid_${IC}_${eos}_errorconvergence.tex"
	set multiplot layout 2,2
	set size square
	set border lw 3
	set tics nomirror
	
	set ylabel 'Density'
	set xrange [0:1]
	set xlabel 'x'
	set key top left
	plot "./../../exact_riemann_solver_idealgas/${IC}.dat" w l ls 2 lw 2 lc 0  title "Exact", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_10_final.dat" u 1:2 w l ls 1 lw 2  title "N = 10", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_40_final.dat" u 1:2 w l ls 1 lw 2 lc 2  title "N = 40", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_160_final.dat" u 1:2 w l ls 1 lw 2 lc 3  title "N = 160", "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_ideal_640_final.dat" u 1:2 w l ls 1 lw 2 lc 4  title "N = 640"



	set logscale x
	set logscale y
	set xlabel 'Number of grid cells'
	set xrange [7:1500]
	set key top right

	set ylabel 'Density L1 error'
	set yrange [0.01:100]
	set ytics ("\$10^{-4}\$"0.0001,"\$10^{-3}\$"0.001, "\$10^{-2}\$" 0.01,"\$10^{-1}\$"0.1, "\$10^{0}\$" 1, "\$10^{1}\$" 10, "\$10^{2}\$" 100, "\$10^{3}\$" 1000, "\$10^{4}\$" 10000, "\$10^{5}\$" 100000)

	plot 20*x**(-1) w l lw 3  title '\footnotesize First order convergence', "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 2  title '\footnotesize Godunov', "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_${eos}_collateddensityerrors.dat" w p pt 7 ps 2  title '\footnotesize MUSCL-Hancock'

	set ylabel 'Velocity L1 error'
	set yrange [0.001:10]

	plot 10*x**(-1) w l lw 3  title '\footnotesize First order convergence', "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collatedvelocityerrors.dat" w p pt 7 ps 2  title '\footnotesize Godunov', "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_${eos}_collatedvelocityerrors.dat" w p pt 7 ps 2  title '\footnotesize MUSCL-Hancock'

	set ylabel 'Pressure L1 error'
	set yrange [0.1:1000]

	plot 1000*x**(-1) w l lw 3  title '\footnotesize First order convergence', "./../output/onefluid_Godunov_Exact_idealgas_${IC}_${eos}_collatedpressureerrors.dat" w p pt 7 ps 2  title '\footnotesize Godunov', "./../output/onefluid_MUSCL_Exact_idealgas_${IC}_${eos}_collatedpressureerrors.dat" w p pt 7 ps 2  title '\footnotesize MUSCL-Hancock'

EOF
