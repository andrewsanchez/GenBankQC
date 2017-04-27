#!/bin/bash

#SBATCH --job-name=genbankfilter
#SBATCH --time=24:00:00
#SBATCH --mem=4000
#SBATCH --workdir=/home/aas229/projects/Genbank-Filter
#SBATCH --output=/common/contrib/databases/genbank_04-12/.info/filtering/genbankfilter.out

module load anaconda/3.latest
srun python run.py /common/contrib/databases/genbank_04-12/ --no_update
