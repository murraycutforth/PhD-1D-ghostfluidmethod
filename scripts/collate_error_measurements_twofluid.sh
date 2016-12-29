# Collate all error measurements for a given RS, FS, TC, eos, into one data file for plotting

RS_pure=$1
RS_mixed=$2
FS=$3
TC=$4
eos1=$5
eos2=$6
GFM=$7

cat ./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${TC}_${eos1}-${eos2}_*_finalerror.dat > ./../output/${GFM}_${FS}_${RS_pure}_${RS_mixed}_${TC}_${eos1}-${eos2}_collatedfinalerrors.dat


