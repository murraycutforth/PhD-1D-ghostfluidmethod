# Collate all error measurements for a given RS, FS, TC, eos, into one data file for plotting

RS_pure=Exact_idealgas
RS_mixed=Exact_idealgas
FS=MUSCL
TC=NE1
GFM=( OriginalGFM M_GFM R_GFM P_GFM)

for j in "${GFM[@]}"
do

	cat ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_*_densityerror.dat > ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_collateddensityerrors.dat

	cat ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_*_velocityerror.dat > ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_collatedvelocityerrors.dat

	cat ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_*_pressureerror.dat > ./../output/${j}_${FS}_${RS_pure}_${RS_mixed}_${TC}_collatedpressureerrors.dat

done


