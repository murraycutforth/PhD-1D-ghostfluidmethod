# Take number of cells, Riemann solver type, and flow solver type as arguments
L=$1
RS=$2
FS=$3
allTCs=( TTC1 TTC2 TTC3 TTC4 TTC5 )

for i in "${allTCs[@]}"
do
	filein=./../output/onefluid_${FS}_${RS}_${i}_ideal_${L}_final.dat
	exactsoln=./../../exact_riemann_solver_idealgas/${i}.dat
	fileout=./../output/onefluid_${FS}_${RS}_${i}_ideal_${L}_final.png
	./plot_onefluid_state_withexactsoln_png.sh $filein $exactsoln $fileout
done
