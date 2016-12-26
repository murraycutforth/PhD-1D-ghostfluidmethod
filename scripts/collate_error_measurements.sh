# Collate all error measurements for a given RS, FS, TC, eos, into one data file for plotting

RS=$1
FS=$2
TC=$3
eos=$4

cat ./../output/onefluid_${FS}_${RS}_${TC}_${eos}_*_finalerror.dat > ./../output/onefluid_${FS}_${RS}_${TC}_${eos}_collatedfinalerrors.dat


