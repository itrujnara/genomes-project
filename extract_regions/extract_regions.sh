#!/usr/bin/env bash
#
# args:
# -f INPUT_FASTA
# -g INPUT_GFF
# -t FEATURE_TYPE
# -c CONTEXT_SIZE


while getopts ":f:g:t:c:" opt; do
  case $opt in
	f) in_fasta="$OPTARG" 
    	;;
    	g) in_gff="$OPTARG"
    	;;
	t) feature="$OPTARG"
	;;
	c) context_size="$OPTARG"
	;;	
    	\?) echo "Invalid option -$OPTARG" >&2
    	exit 1
    	;;
  esac

  case $OPTARG in
    	-*) echo "Option $opt needs a valid argument"
    	exit 1
    	;;
  esac
done

# make fai
zcat $in_fasta > temp1.fa
samtools faidx temp1.fa --fai-idx temp2.fai
rm temp1.fa

# extract relevant annotations from GFF
zcat $in_gff | awk -v feat=$feature -F'\t' 'BEGIN { OFS = FS } !/^#/ && $3==feat;' | gzip > temp3.gff.gz

# add context
bedtools slop -i temp3.gff.gz -b $context_size -g temp2.fai | gzip > temp4.gff.gz

# extract regions to TSV
bedtools getfasta -fi $in_fasta -bed temp4.gff.gz | gzip > $(basename $in_fasta).$feature.gff3.gz
