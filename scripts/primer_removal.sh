#!/bin/bash

#SBATCH -A uppmax2025-2-344
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 4:00:00
#SBATCH -J primer_trimming

module load bioinfo-tools
module load R/4.3.1
module load R_packages/4.3.1
module load cutadapt

# Hitta alla innersta fastq-mappar
fastq_dirs=$(find /home/haddock/private/järv/data \
    -type d -path "*02-FASTQ/250903_VH00203_554_AAH7JV3M5")

for d in $fastq_dirs; do
    echo "Kör trimming på: $d"
    Rscript primer_removal.R "$d"
done
