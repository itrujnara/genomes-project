#!/usr/bin/env bash

for i in $(ls $1)
do
	basename="${i##*_}"
	path=$(realpath $1/$i/${basename}.fasta.gz)
	echo "$i,$path"
done
