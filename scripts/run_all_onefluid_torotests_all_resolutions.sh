# Compile code
cd ..
make
cd ./scripts

# Exit if any simulations fail
set -e

# First set up the settings file (leave RS_pure, FS, CFL unchanged)
sed -i 's/numGC.*/numGC 2/g' ./../settings_file.txt
sed -i 's/fluid1_gamma.*/fluid1_gamma 1.4/g' ./../settings_file.txt
sed -i 's/eos1.*/eos1 ideal/g' ./../settings_file.txt
sed -i 's/BC_L.*/BC_L transmissive/g' ./../settings_file.txt
sed -i 's/BC_R.*/BC_R transmissive/g' ./../settings_file.txt
sed -i 's/sim.*/sim onefluid/g' ./../settings_file.txt

# List of all the test cases
alltcs=( TTC1 TTC2 TTC3 TTC4 TTC5 )

# List of all resolutions to be used
allLs=( 10 20 40 80 160 320 640 1280 )

for i in "${alltcs[@]}"
do
	sed -i "s/IC.*/IC ${i}/g" ./../settings_file.txt
	
	for j in "${allLs[@]}"
	do
		sed -i "s/length.*/length ${j}/g" ./../settings_file.txt
		cd ..
		./1D_Euler_GFM.exe
		cd ./scripts
	done
done
