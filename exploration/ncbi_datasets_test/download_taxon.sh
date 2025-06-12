#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "Usage: $0 <taxon_id>"
    exit 1
fi

taxon_id=$1
datasets summary genome taxon --assembly-level complete,chromosome --reference --annotated $taxon_id > ${taxon_id}.json
python3 get_genome_ids.py ${taxon_id}.json > ${taxon_id}.txt
shuf ${taxon_id}.txt | head -n 100 > ${taxon_id}_100.txt
datasets download genome accession --inputfile ${taxon_id}_100.txt --filename ${taxon_id}_100.zip --include genome,gff3 --dehydrated
unzip ${taxon_id}_100.zip -d ${taxon_id}_100
