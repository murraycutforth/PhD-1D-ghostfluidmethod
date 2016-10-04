filein=$1
fileout=$2

gnuplot <<- EOF
	set terminal postscript eps enhanced color font 'Palatino,60' size 15,10

	set output '$fileout'
	set border lw 3
	set xlabel 'Time'
	set ylabel 'Conservation error'
	
	plot "$filein" u 1:2 w lp ps 4 lw 3 pt 4 title 'Mass', "$filein" u 1:3 w lp ps 4 lw 3 pt 4 title 'Momentum', "$filein" u 1:4 w lp ps 4 lw 3 pt 4 title 'Energy'

EOF
