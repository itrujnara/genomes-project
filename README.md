# Genomes ML project

This repository contains files used in the LLM genome annotation project and the related manuscript.

# Directories and files

## Directories used for the manuscript
The following directories contained in this repository are directly relevant to the manuscript:
- `infer_chromosome` – files related to the chromosome 21 inference run and its analysis
- `notebooks` – Jupyter notebooks used for data postprocessing and plotting
- `scripts` – basic scripts for Evo and Evo 2 inference, used for the PoC gene annotation analysis

## Directories related to the manuscript
The following directories are connected to the manuscript, but not directly relevant:
- `extract_regions` – a test of region extraction from FASTA based on GFF
- `gene_boundaries` – a test of gene boundary inference with Evo
- `minihyena` – a test of a small model based on Evo
- `test_geneid` – a test of using Evo output to improve geneid

## Other directories
The following other directories are not related to the manuscript or purely technical:
- `exploration` – exploratory analyses done before the Evo analysis; not fully maintained
- `logs` – SLURM job stdout and stderr, not tracked

Most of the directories contain their own README with details.
