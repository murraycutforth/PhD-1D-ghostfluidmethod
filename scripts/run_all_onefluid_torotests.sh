# Compile code
cd ..
make
cd ./scripts

# Exit if any simulations fail
set -e

# First set up the settings file (leave RS_pure, FS, CFL unchanged)
sed -i 's/length.*/length 100/g' ./../settings_file.txt
sed -i 's/numGC.*/numGC 2/g' ./../settings_file.txt
sed -i 's/fluid1_gamma.*/fluid1_gamma 1.4/g' ./../settings_file.txt
sed -i 's/eos1.*/eos1 ideal/g' ./../settings_file.txt
sed -i 's/BC_L.*/BC_L transmissive/g' ./../settings_file.txt
sed -i 's/BC_R.*/BC_R transmissive/g' ./../settings_file.txt
sed -i 's/sim.*/sim onefluid/g' ./../settings_file.txt

# Now create list of all the test cases and run for each one
alltcs=( TTC1 TTC2 TTC3 TTC4 TTC5 )

for i in "${alltcs[@]}"
do
	sed -i "s/IC.*/IC ${i}/g" ./../settings_file.txt
	cd ..
	./1D_Euler_GFM.exe
	cd ./scripts
done
