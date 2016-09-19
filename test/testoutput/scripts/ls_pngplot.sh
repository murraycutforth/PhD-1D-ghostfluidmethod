filein=$1
fileout=$2

gnuplot <<- EOF
	set terminal png font 'arial,30' size 1500,1000
	set output '$fileout'
	set border lw 3
	unset key
	set xlabel 'x'
	set ylabel 'level set function'
	set yrange [-1:1]
	plot '$filein' u 1:2 w lp ps 4 lw 3 pt 4
EOF
