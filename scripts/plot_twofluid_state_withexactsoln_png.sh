step=final
N=100
GFM=OriginalGFM
FS=MUSCL
RS_pure=Exact_idealgas
RS_mixed=Exact_idealgas
IC=NE1
fileout=${GFM}_${IC}_realfluidstate.tex
filein1=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${N}_realfluid1_${step}.dat
filein2=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${N}_realfluid2_${step}.dat
filein3=./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${IC}_${N}_cellwiseerror.dat
exactsoln=$1

gnuplot <<- EOF
	set terminal epslatex color size 16cm,20cm
	set output '$fileout'
	set border lw 3
	unset key
	set xtics nomirror
	set multiplot layout 4,2 noenhanced
	
	set lmargin at screen 0.1; set rmargin at screen 0.45; set tmargin at screen 0.99; set bmargin at screen 0.69;
	set ylabel 'density'
	unset xtics
	plot "$exactsoln" u 1:2 w l lw 3 lc 4, "$filein1" u 1:2 w p ps 2 pt 4 lw 3 lc 1, "$filein2" u 1:2 w p ps 2 pt 4 lw 3 lc 3

	set ylabel 'velocity'
	set lmargin at screen 0.64; set rmargin at screen 0.99; set tmargin at screen 0.99; set bmargin at screen 0.69;
	plot "$exactsoln" u 1:3 w l lw 3 lc 4, "$filein1" u 1:3 w p ps 2 pt 4 lw 3 lc 1, "$filein2" u 1:3 w p ps 2 pt 4 lw 3 lc 3

	set ylabel 'density error'
	set xlabel 'x'
	set xtics
	set xtics nomirror
	set lmargin at screen 0.1; set rmargin at screen 0.45; set tmargin at screen 0.68; set bmargin at screen 0.58;
	plot "$filein3" u 1:2 w filledcurves

	set ylabel 'velocity error'
	set lmargin at screen 0.64; set rmargin at screen 0.99; set tmargin at screen 0.68; set bmargin at screen 0.58;
	plot "$filein3" u 1:3 w filledcurves

	set ylabel 'pressure'
	unset xtics
	unset xlabel
	set lmargin at screen 0.1; set rmargin at screen 0.45; set tmargin at screen 0.51; set bmargin at screen 0.21;
	plot "$exactsoln" u 1:4 w l lw 3 lc 4, "$filein1" u 1:5 w p ps 2 pt 4 lw 3 lc 1, "$filein2" u 1:5 w p ps 2 pt 4 lw 3 lc 3

	set ylabel 'specific internal energy'
	set lmargin at screen 0.64; set rmargin at screen 0.99; set tmargin at screen 0.51; set bmargin at screen 0.21;
	plot "$exactsoln" u 1:5 w l lw 3 lc 4, "$filein1" u 1:6 w p ps 2 pt 4 lw 3 lc 1, "$filein2" u 1:6 w p ps 2 pt 4 lw 3 lc 3
	
	set ylabel 'pressure error'
	set xlabel 'x'
	set xtics nomirror
	set lmargin at screen 0.1; set rmargin at screen 0.45; set tmargin at screen 0.2; set bmargin at screen 0.1;
	plot "$filein3" u 1:4 w filledcurves

	set ylabel 'specific internal energy error'
	set lmargin at screen 0.64; set rmargin at screen 0.99; set tmargin at screen 0.2; set bmargin at screen 0.1;
	plot "$filein3" u 1:5 w filledcurves
EOF
