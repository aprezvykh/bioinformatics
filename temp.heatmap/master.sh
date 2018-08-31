#!/bin/bash
seq_num=1000
echo "Total cycles is set to $seq_num"
for f in `seq $seq_num`;
do
			sensors | grep "Core" > temp.txt
		sensors | grep "Core" > $(date -d "today" +"%Y%m%d%H%M%S").txt
		Rscript hm.R
		rm temp.txt
	sleep 2
	echo "Cycle $f completed"
	num=$(($seq_num - $f))
	echo "$num cycles is left"

done
echo "Convertation started!"
convert -delay 2 -loop 0 *.png heat.gif
rm *.png
Rscript plot.R
fld=$(date -d "today" +"%Y%m%d%H%M")
mkdir $fld
echo "Removing some..."
cp *.png *.gif $fld/
rm *.png *.gif *.txt
echo "done!"
