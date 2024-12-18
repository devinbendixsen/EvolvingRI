#!/bin/bash

while read bam vcf; do
  echo "$bam"
  sbatch run_freebayes.sh $bam $vcf
done < Hybrid_Dynamics_freebayes_samples.txt
