basefilename=$1
i=0
increment=1

# Create png plot of each output
while [ -e ./../output/${basefilename}${i}.dat ]; do
	filein=./../output/${basefilename}${i}.dat
	fileout=./../output/videoframes/${basefilename}${i}.png
	./state_pngplot.sh $filein $fileout
	let i=i+increment
done

# Create video of output
avconv -r 25 -i ./../output/videoframes/${basefilename}%d.png ./../output/processed/${basefilename}.mp4

# Open video in vlc
if [ -e ./../output/processed/${basefilename}.mp4 ]
then
	vlc ./../output/processed/${basefilename}.mp4
fi
