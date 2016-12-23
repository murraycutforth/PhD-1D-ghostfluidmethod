filein=$1
exactsoln=$2
fileout=$3

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,30' size 3000,2000
	set output '$fileout'
	set border lw 5
	unset key
	set multiplot layout 2,2 title "$filein" noenhanced
	set xrange [0:1]
	
	set xlabel 'x'
	set ylabel 'density'
	plot "$exactsoln" u 1:2 w l lw 5 lc 3, "$filein" u 1:2 w p ps 3 lw 3 pt 4 lc 1

	set ylabel 'velocity'
	plot "$exactsoln" u 1:3 w l lw 5 lc 3, "$filein" u 1:3 w p ps 3 lw 3 pt 4 lc 1

	set ylabel 'pressure'
	plot "$exactsoln" u 1:4 w l lw 5 lc 3, "$filein" u 1:5 w p ps 3 lw 3 pt 4 lc 1

	set ylabel 'specific internal energy'
	plot "$exactsoln" u 1:5 w l lw 5 lc 3, "$filein" u 1:6 w p ps 3 lw 3 pt 4 lc 1
EOF
