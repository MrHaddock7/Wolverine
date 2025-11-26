#!/bin/bash
#SBATCH -A uppmax2025-2-344
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 10:00:00
#SBATCH -J picrust_job

module load Miniconda3

conda activate picrust2

picrust2_pipeline.py -s /home/haddock/private/Wolverine/data/ASV_abundance_sequences/rep_seqs.fna -i /home/haddock/private/Wolverine/data/ASV_abundance_sequences/table.biom -o /home/haddock/private/Wolverine/results/picrust2_output_26_11_2025/picrust2_out_pipeline_run -p 20