#!/usr/bin/env bash

#SBATCH --output=/users/cn/itrujnara/genomes-project/logs/%j.out
#SBATCH --error=/users/cn/itrujnara/genomes-project/logs/%j.err
#
#SBATCH --time=10
#SBATCH --mem=24G
#SBATCH --qos=normal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:3g.40gb

# basic setup
set -euo pipefail

# environment
module load Python
module load CUDA
source /users/cn/itrujnara/envs/venv_gpu/bin/activate

# script
python3 evo_inference.py $*
