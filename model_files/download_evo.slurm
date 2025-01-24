#!/bin/bash

####################
# slurm directives #
####################

# where to put stdout / stderr
#SBATCH --output=./logs/%j.out
#SBATCH --error=./logs/%j.err

# time limit in minutes
#SBATCH --time=5

# which qos to use (similar to queue in Grid Engine)
#SBATCH --qos=vshort
#SBATCH --partition=gpu
#SBATCH --gres=gpu:3g.40gb:1

# memory (MB)
#SBATCH --mem=32G

###################
# set environment #
###################
module load Python
module load CUDA
source /users/cn/itrujnara/venv_gpu/bin/activate

###############
# run command #
###############
python3 download_evo.py
