#!/bin/bash

# Hitta alla innersta fastq-mappar
fastq_dirs=$(find /Users/william/Library/CloudStorage/OneDrive-Uppsalauniversitet/Slutkurs/OneDrive_1_2025-11-11 \
    -type d -path "*02-FASTQ/250903_VH00203_554_AAH7JV3M5")

for d in $fastq_dirs; do
    echo "Kör trimming på: $d"
    Rscript primer_removal.r "$d"
done
