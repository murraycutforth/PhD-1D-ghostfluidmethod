filein=$1
fileout=$2
exact=$3

gnuplot <<- EOF
	set terminal postscript eps enhanced color font 'Palatino,60' size 30,20

	set output '$fileout'
	set border lw 3
	unset key
	set multiplot layout 2,2
	set xrange [0:1]
	
	set xlabel 'x'
	set ylabel 'density'
	plot "$exact" u 1:2 w l lc 3 lw 6, "$filein" u 1:2 w lp ps 4 lw 3 pt 4 lc 1

	set ylabel 'velocity'
	plot "$exact" u 1:3 w l lc 3 lw 6, "$filein" u 1:3 w lp ps 4 lw 3 pt 4 lc 1

	set ylabel 'pressure'
	plot "$exact" u 1:4 w l lc 3 lw 6, "$filein" u 1:5 w lp ps 4 lw 3 pt 4 lc 1

	set ylabel 'specific internal energy'
	plot "$exact" u 1:5 w l lc 3 lw 6, "$filein" u 1:6 w lp ps 4 lw 3 pt 4 lc 1
EOF
