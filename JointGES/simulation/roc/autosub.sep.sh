#!/usr/bin/env bash

for file in `ls ../data/*neigh_2*.rda`
do
	#for lambda in 0.125 0.25 0.5 1 3 5 7 9 11 15
	for lambda in 1.5 2 4
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${lambda}" < sep.R" > tmp.sh
		cat template.sep.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
