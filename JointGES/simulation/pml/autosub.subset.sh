#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for lambda in 0.3 1 3 10
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${lambda}" < subset.R" > tmp.sh
		cat template.joint.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
