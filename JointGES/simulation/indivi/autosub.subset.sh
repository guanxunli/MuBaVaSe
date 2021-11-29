#!/usr/bin/env bash

#for file in `ls ../data/*.rda`
for file in `ls ../data/p_100_neigh_2_n_900_k_5.rda`
do
	#for lambda in 2 3 4 5
	#for lambda in 1
	for lambda in 2
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${lambda}" < subset.R" > tmp.sh
		cat template.joint.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
