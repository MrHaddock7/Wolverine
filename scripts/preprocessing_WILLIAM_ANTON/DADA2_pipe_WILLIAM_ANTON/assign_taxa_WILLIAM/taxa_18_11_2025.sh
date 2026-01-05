#!/bin/bash -l
#SBATCH -A uppmax2025-2-344       # Replace with your project, e.g. sens2024012
#SBATCH -p core                  # Use the 'core' partition for normal CPU jobs
#SBATCH -n 4                     # Number of cores
#SBATCH -t 2:00:00              # Wall time (2 hours here)
#SBATCH -J taxa_analysis        # Job name
#SBATCH -o taxa_%j.out      # STDOUT (%j = job ID)
#SBATCH -e taxa_%j.err      # STDERR

# --- Load R environment ---

module load R/4.3.1
module load R_packages/4.3.1

Rscript /home/haddock/private/Wolverine/scripts/assign_taxa/assign_taxa.R 
