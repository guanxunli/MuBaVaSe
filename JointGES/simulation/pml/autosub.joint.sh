#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for lambda in 1 2 3 4 5
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${lambda}" < joint.R" > tmp.sh
		cat template.joint.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
