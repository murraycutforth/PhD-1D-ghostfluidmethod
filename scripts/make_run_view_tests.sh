basefilename=$1
i=0
increment=1

# Delete old data files
rm ../test/testoutput/*.dat
rm ../test/testoutput/videoframes/*.png

# Make and run test code
./../make
./../make tester
./../test/maintest

# Create png plot of each output
while [ -e ./../test/testoutput/${basefilename}${i}.dat ]; do
	filein=./../test/testoutput/${basefilename}${i}.dat
	fileout=./../test/testoutput/videoframes/${basefilename}${i}.png
	./state_pngplot.sh $filein $fileout
	let i=i+increment
done

# Create video of output
avconv -r 25 -i ./../test/testoutput/videoframes/${basefilename}%d.png ./../test/testoutput/videoframes/${basefilename}.mp4

# Open video in vlc
vlc ./../test/testoutput/videoframes/${basefilename}.mp4

