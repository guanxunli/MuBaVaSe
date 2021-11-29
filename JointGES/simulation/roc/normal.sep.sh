#!/usr/bin/env bash

for file in `ls ../data/*neigh_2*.rda`
do
	#for lambda in 0.125 0.25 0.5 1 3 5 7 9 11 15
	for lambda in 1.5
	do
		file=`basename $file`
		/home/yuhaow/software/R/bin/R --no-save --args ../data/${file} ${lambda} < sep.R
	done
done
