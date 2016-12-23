filein=$1
fileout=$2

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,30' size 3000,2000
	set output '$fileout'
	set border lw 5
	unset key
	set multiplot layout 2,2 title "$filein" noenhanced
	
	set xlabel 'x'
	set ylabel 'density'
	plot "$filein" u 1:2 w lp ps 4 lw 3 pt 4

	set ylabel 'velocity'
	plot "$filein" u 1:3 w lp ps 4 lw 3 pt 4

	set ylabel 'pressure'
	plot "$filein" u 1:5 w lp ps 4 lw 3 pt 4

	set ylabel 'specific internal energy'
	plot "$filein" u 1:6 w lp ps 4 lw 3 pt 4
EOF
