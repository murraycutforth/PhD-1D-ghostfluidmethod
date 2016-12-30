step=$1
N=$3
GFM=R_GFM
FS=MUSCL
RS_pure=HLLC_idealgas
RS_mixed=M_HLLC
IC=TTC1
eos=ideal-ideal
fileout=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluidstate.png
filein1=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluid1_${step}.dat
filein2=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${eos}_${N}_realfluid2_${step}.dat
exactsoln=$2

gnuplot <<- EOF
	set terminal png enhanced font 'Palatino,30' size 3000,2000
	set output '$fileout'
	set border lw 5
	unset key
	set multiplot layout 2,2 title "$fileout" noenhanced
	
	set xlabel 'x'
	set ylabel 'density'
	plot "$exactsoln" u 1:2 w l lw 3 lc 4, "$filein1" u 1:2 w p ps 4 pt 4 lw 3 lc 1, "$filein2" u 1:2 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'velocity'
	plot "$exactsoln" u 1:3 w l lw 3 lc 4, "$filein1" u 1:3 w p ps 4 pt 4 lw 3 lc 1, "$filein2" u 1:3 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'pressure'
	plot "$exactsoln" u 1:4 w l lw 3 lc 4, "$filein1" u 1:5 w p ps 4 pt 4 lw 3 lc 1, "$filein2" u 1:5 w p ps 4 pt 4 lw 3 lc 3

	set ylabel 'specific internal energy'
	plot "$exactsoln" u 1:5 w l lw 3 lc 4, "$filein1" u 1:6 w p ps 4 pt 4 lw 3 lc 1, "$filein2" u 1:6 w p ps 4 pt 4 lw 3 lc 3
EOF
