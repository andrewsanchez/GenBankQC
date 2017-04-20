#!/bin/bash

#SBATCH --job-name=ncbitk
#SBATCH --output=/scratch/aas229/ncbitk.out
#SBATCH --time=72:00:00
#SBATCH --workdir=/home/aas229/projects/NCBITK

module load anaconda/3.latest

srun python run.py /common/contrib/databases/genbank_04-12 --update
