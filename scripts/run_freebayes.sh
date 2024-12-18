#!/bin/bash

#Submit this script with: sbatch thefilenamex

#SBATCH -A naiss2023-22-62
#SBATCH --time=12:00:00   # walltime
#SBATCH -n 2
#SBATCH -J "freebayes_Hybrid"   # job name
#SBATCH --mail-user=devin.bendixsen@zoologi.su.se   # email address

# conda activate freebayes

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

freebayes -f ../data/reference_genome/genome.fa ../results/bam_files/$1 >../results/vcfs/$2
