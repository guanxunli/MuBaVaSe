#!/usr/bin/env bash

for file in `ls ../data/*neigh_2*.rda`
do
	#for lambda in 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5
	for lambda in 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} 3 ${lambda}" < subset.R" > tmp.sh
		cat template.joint.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
