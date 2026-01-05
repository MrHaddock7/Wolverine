#!/bin/bash

#SBATCH -A uppmax2025-2-344
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 14:00:00
#SBATCH -J primer_trimming
#SBATCH -o primer_%j.out      # STDOUT (%j = job ID)
#SBATCH -e primer_%j.err      # STDERR

module load bioinfo-tools
module load R/4.3.1
module load R_packages/4.3.1
module load cutadapt

# Hitta alla innersta fastq-mappar
fastq_dirs=$(find /proj/rnaseq01/viroids/Viroid-transcriptome-mining/slutkurs/G.gulo_16S_zoo \
    -type d -path "*02-FASTQ/250903_VH00203_554_AAH7JV3M5")

for d in $fastq_dirs; do
    echo "Kör trimming på: $d"
    Rscript /home/haddock/private/järv/scripts/primer_removal.R "$d"
done
