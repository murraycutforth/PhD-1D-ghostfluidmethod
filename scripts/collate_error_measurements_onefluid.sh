# Collate all error measurements for a given RS, FS, TC, eos, into one data file for plotting

RS=Exact_idealgas
FS=MUSCL
TC=( TTC1 TTC2 TTC3 TTC4 TTC5 GDA )

for j in "${TC[@]}"
do
	cat ./../output/onefluid_${FS}_${RS}_${j}_ideal_*_densityerror.dat > ./../output/onefluid_${FS}_${RS}_${j}_ideal_collateddensityerrors.dat
	cat ./../output/onefluid_${FS}_${RS}_${j}_ideal_*_velocityerror.dat > ./../output/onefluid_${FS}_${RS}_${j}_ideal_collatedvelocityerrors.dat
	cat ./../output/onefluid_${FS}_${RS}_${j}_ideal_*_pressureerror.dat > ./../output/onefluid_${FS}_${RS}_${j}_ideal_collatedpressureerrors.dat
done


