# Pathway Heatmap Guide

This script (`scripts/function_analysis/pathway_heatmap_top20.R`) builds a ggpicrust2 heatmap for the top 20 MetaCyc pathways and flags which are significantly different across zoos.

## How the code works
- Load metadata (`data/metadata.csv`) and PICRUSt2 pathway abundances (`results/picrust2_output_26_11_2025/.../path_abun_unstrat.tsv.gz`), then keep only overlapping samples so counts and metadata align.
- Pick the 20 most abundant MetaCyc pathways across all samples; this focuses the heatmap on dominant functions.
- Run `pathway_daa(..., daa_method = "LinDA", group = "Zoo")` to test pathway differences between zoos. The script adapts to ggpicrust2 returning either a data.frame or a list and uses `p_adjust` (BH q-values) to flag significant pathways.
- Add `[sig]` to any pathway with `p_adjust < 0.05`, plot with `pathway_heatmap()`, and save two outputs: `results/top20_pathway_heatmap.png` and the differential table `results/top20_pathway_daa.csv`.
- Run PCA on the top-20 matrix (samples as rows) to see whether functional profiles separate by `Zoo`; save as `results/top20_pathway_pca.png` with variance explained on axes.
- Run PERMANOVA (Bray-Curtis on the top-20 matrix) to test om Zoo förklarar variansen; results till `results/top20_pathway_permanova.txt`.

## How to run
```sh
Rscript scripts/function_analysis/pathway_heatmap_top20.R
```
Update the file paths or group column in the script if your metadata or PICRUSt outputs live elsewhere.

## What the heatmap shows
Rows are pathways, columns are samples; darker fill indicates higher abundance. `[sig]` marks pathways whose abundance shifts across zoos (LinDA q < 0.05), highlighting the “useful” functional signals to investigate first.

## What the PCA shows
Each point is a sample; proximity means similar KO-based pathway profiles within the top-20 set. Separation by color (Zoo) suggests site-specific functional signatures; dashed ellipses encircle each zoo to make group differences easier to see. Axis labels show variance explained by PC1/PC2.

## What PERMANOVA shows
Statistisk test av om Zoo förklarar variationen i pathway-profilerna (Bray-Curtis). F- och p-/q-värden ger stöd för om skillnaderna är signifikanta; ligger i `results/top20_pathway_permanova.txt`.

## Use case for wolverine microbiota
The plot highlights functional capacities that differ between zoo populations—helpful for linking diet, habitat, or health status to metabolic potential. Track whether husbandry changes move key pathways toward desired profiles or reveal site-specific dysbiosis worth addressing.
