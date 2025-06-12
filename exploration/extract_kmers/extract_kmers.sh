#!/usr/bin/env bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_fasta> <kmer_size>"
    exit 1
fi

input_fasta=$1
base_filename=$(basename "$input_fasta" | sed 's/\.[^.]*\.[^.]*$//' | sed 's/_genomic//')
jellyfish count -m $2 -s 100M -t 10 -o "${base_filename}_$2mers.jf" <(zcat "$input_fasta")
echo -e "kmer\t${base_filename}" > "${base_filename}_$2mers.tsv"
jellyfish dump "${base_filename}_$2mers.jf" | tr '\n' '\t' | tr '>' '\n' | sed 's/\t$//g' | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $2,$1}' | tail -n +2 >> "${base_filename}_$2mers.tsv"
