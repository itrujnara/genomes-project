#!/usr/bin/env bash

sbatch infer_chunk.sbatch --parquet ch38_chr21_chunks_1mbp.pq --id chunk$1 --prefix ch38_chr21
