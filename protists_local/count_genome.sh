#!/usr/bin/env bash

filename="$1"
basename="${filename%.*}"

zcat $1 > tmp
jellyfish count -m 5 -s 1M -t 12 -o "${basename}.jf" tmp
rm tmp
