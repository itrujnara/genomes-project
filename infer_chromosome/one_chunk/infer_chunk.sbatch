#!/usr/bin/env bash

#SBATCH --output=/users/cn/itrujnara/genomes-project/logs/%j.out
#SBATCH --error=/users/cn/itrujnara/genomes-project/logs/%j.err
#
#SBATCH --time=30
#SBATCH --mem=24G
#SBATCH --qos=normal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:3g.40gb

# basic setup
set -euo pipefail

# environment
module load Python/3.11.5
module load CUDA
module load cuDNN
source /users/cn/itrujnara/envs/evo2/bin/activate

# env variables
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# script
python3 /users/cn/itrujnara/genomes-project/infer_chromosome/one_chunk/infer_chunk.py $*
