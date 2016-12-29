step=$1
N=160
GFM=OriginalGFM
FS=Godunov
RS_pure=Exact_idealgas
RS_mixed=Exact_idealgas
IC=TTC1
eos=ideal-ideal
fileout=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluidstate.png
filein1=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluid1_${step}.dat
filein2=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluid2_${step}.dat

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,30' size 3000,2000
	set output '$fileout'
	set border lw 5
	unset key
	set multiplot layout 2,2 title "$fileout" noenhanced
	
	set xlabel 'x'
	set ylabel 'density'
	plot "$filein1" u 1:2 w p ps 4 pt 4 lw 3, "$filein2" u 1:2 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'velocity'
	plot "$filein1" u 1:3 w p ps 4 pt 4 lw 3, "$filein2" u 1:3 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'pressure'
	plot "$filein1" u 1:5 w p ps 4 pt 4 lw 3, "$filein2" u 1:5 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'specific internal energy'
	plot "$filein1" u 1:6 w p ps 4 pt 4 lw 3, "$filein2" u 1:6 w p ps 4 pt 4 lw 3 lc 3
EOF
