#!/bin/bash
#SBATCH -A uppmax2025-2-344
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 10:00:00
#SBATCH -J picrust_job

module load Miniconda3

echo "Loading conda"

source /sw/arch/eb/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate picrust2

echo "Starting picrust"

picrust2_pipeline.py -s /home/haddock/private/Wolverine/data/ASV_abundance_sequences/rep_seqs.fna -i /home/haddock/private/Wolverine/data/ASV_abundance_sequences/table.biom -o /home/haddock/private/Wolverine/results/picrust2_output_26_11_2025/picrust2_out_pipeline_run -p 20