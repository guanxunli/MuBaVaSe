#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for lambda in 0.5 1 2 3 4 5
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${lambda}" < sep.R" > tmp.sh
		cat template.sep.sh tmp.sh > scripts/${file}.${lambda}.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sh
	done
done
