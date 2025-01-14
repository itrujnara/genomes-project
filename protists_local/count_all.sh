#!/usr/bin/env bash

for i in $(ls $1)
do
	basename=${i##*_}
	fastaname="$1/$i/${basename}.fasta.gz"
	zcat $fastaname > "$1/$i/tmp"
	jellyfish count -m 5 -s 1M -t 12 -o "$1/$i/${basename}.jf" $1/$i/tmp
	rm $1/$i/tmp
done
