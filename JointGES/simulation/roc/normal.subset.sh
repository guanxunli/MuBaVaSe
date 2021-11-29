#!/usr/bin/env bash

for file in `ls ../data/*neigh_2*.rda`
do
	#for lambda in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
	for lambda in 0.01 0.05
	do
		file=`basename $file`
		/home/yuhaow/software/R/bin/R --no-save --args ../data/${file} 2 ${lambda} < subset.R
	done
done
