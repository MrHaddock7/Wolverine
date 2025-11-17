#!/bin/bash -l
#SBATCH -A uppmax2025-2-344       # Replace with your project, e.g. sens2024012
#SBATCH -p core                  # Use the 'core' partition for normal CPU jobs
#SBATCH -n 10                     # Number of cores
#SBATCH -t 24:00:00              # Wall time (2 hours here)
#SBATCH -J analysis_gm_w        # Job name
#SBATCH -o logs/taxa_%j.out      # STDOUT (%j = job ID)
#SBATCH -e logs/taxa_%j.err      # STDERR

# --- Load R environment ---

module load R/4.3.1
module load R_packages/4.3.1

# -- Create a run directory where all checkpoints are saved
RUN_DIR="/home/haddock/private/Wolverine/runs/run_trimming_taxa_17112025_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR"
cp test_taxa.R savepoint1.RData "$RUN_DIR"/
cd "$RUN_DIR"

mkdir -p logs

# -- Run the scripts: 
Rscript /home/haddock/private/Wolverine/scripts/trimming/filtering.R 

Rscript /home/haddock/private/Wolverine/scripts/trimming/error_and_chimeras.R 

Rscript /home/haddock/private/Wolverine/scripts/assign_taxa/assign_taxa.R 
