#!/usr/bin/env bash
datasets summary gene taxon 9606 --as-json-lines | dataformat tsv gene --fields gene-id,gene-type --elide-header > human_genes.txt
