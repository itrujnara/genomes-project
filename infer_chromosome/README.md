# Infer chromosome

This directory contains the code and part of the data for the human chromosome 21 analysis.

## Directories
This directory contains the following subdirectories:
- `benchmark` – code for benchmarking Evo 2 inference efficiency with different options
- `ch38_data` – data for the human genome obtained from NCBI; not tracked
- `compare_gff` – annotation data for the feature IC analysis in GFF and database format
- `conservation` – data for the comparison between Evo and phastCons
- `one_chunk` – data and code for the inference parameter analyses
- `run_pipeline` – tests of a Nextflow pipeline for Evo inference; the pipeline could not be finished
- `whole_chromosome` – code for running Evo 2 inference on chromosome 21

## Files
This directory contains the following files:
- `test.py` – a prototype script to test Evo 2 inference
- `test.sbatch` – a SLURM submit script for the Python script

## Notes
Some data files (e.g. human genome FASTAs) are not tracked due to their size.
In a reproduction attempt, the correct files must be provided.

The directory contains multiple SLURM submit scripts. 
They will likely need modification to work on any specific SLURM cluster.
