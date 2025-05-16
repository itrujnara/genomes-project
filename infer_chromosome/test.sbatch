#!/usr/bin/env bash
#
#SBATCH --output /users/cn/itrujnara/genomes-project/logs/%j.out
#SBATCH --error /users/cn/itrujnara/genomes-project/logs/%j.err
#SBATCH --partition gpu
#SBATCH --gres gpu:3g.40gb
#SBATCH --mem 32G
#SBATCH --time 10

module load Python/3.11.5
module load CUDA
module load cuDNN

source /users/cn/itrujnara/envs/evo2/bin/activate

python3 /users/cn/itrujnara/genomes-project/infer_chromosome/test.py
