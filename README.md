# Genomes ML project

This repository contains files used locally in the eukaryotic genomes ML project. 

It contains the following directories and files:
- `explore_kmers` – an R Shiny application for exploring k-mer spectra generated with Jellyfish
- `extract_kmers` – scripts and sample files for testing k-mer spectrum generation with Jellyfish
- `extract_regions` – testing sequence extraction from a genome based on a GFF wiht AGAT
- `get_kmers` – get kmers in genomes with an end-to-end pipeline
- `hyena` – a Hyena config
- `logs` – SLURM job outputs, not tracked
- `model_files` – initial test of loading the Evo model
- `ncbi_datasets_test` – scripts to test semi-automated data download from NCBI with the new Datasets CLI; the downloaded data is excluded from version control due to excessive size
- `notebooks` – Jupyter notebooks to test models and analyze outputs
- `protists_local` – some FASTA, GFF and BED files of protist genomes to test certain k-mer commands
- `run_stimulus` – test Stimulus tools on the minihyena model
- `test_geneid` – test geneid with LLM hints
- `mount.sh` – a script to mount genome directories from the CRG cluster inside this working directory
