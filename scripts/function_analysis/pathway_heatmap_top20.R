#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(vegan)
  library(tidyr)
})
source("config.R")

metadata_path <- file.path(home, "data", "metadata.csv")
pathway_abun_path <- file.path(
  home,
  "results/picrust2_output_26_11_2025/picrust2_out_pipeline_run/pathways_out/path_abun_unstrat.tsv.gz"
)

metadata <- read.csv(metadata_path, sep = ";", fileEncoding = "latin1")
if (!"NGI.ID" %in% names(metadata)) {
  stop("metadata.csv must contain an NGI.ID column")
}
rownames(metadata) <- metadata$NGI.ID

metacyc_abundance <- readr::read_tsv(pathway_abun_path, show_col_types = FALSE)
if (!"pathway" %in% names(metacyc_abundance)) {
  stop("path_abun_unstrat.tsv.gz must contain a 'pathway' column")
}
pathway_abundance <- metacyc_abundance %>% tibble::column_to_rownames("pathway")


# Låg statistisk kraft för vissa pathways som finns i < 3 samples (5 för tuffare filtrering)

keep <- rowSums(pathway_abundance > 0) >=3
pathway_abundance <- pathway_abundance[keep, , drop = FALSE]


# Metadata table should only include the samples from the picrust output

common_samples <- intersect(colnames(pathway_abundance), rownames(metadata))

pathway_abundance <- pathway_abundance[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Differential abundance to flag pathways that differ by Zoo
daa_res <- pathway_daa(
  abundance = pathway_abundance,
  metadata = metadata,
  group = "Zoo",
  daa_method = "LinDA"
)

# PCoA on filtered abundance (Bray-Curtis), colored by Zoo
dist_bray <- vegan::vegdist(t(pathway_abundance), method = "bray")
pcoa_bray <- stats::cmdscale(dist_bray, k = 2, eig = TRUE)

pcoa_bray_df <- as.data.frame(pcoa_bray$points) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata %>% rownames_to_column("SampleID"), by = "SampleID")

eig_pct_bray <- round(100 * pcoa_bray$eig / sum(pcoa_bray$eig[pcoa_bray$eig > 0]), 1)

pcoa_bray_plot <- ggplot(pcoa_bray_df, aes(x = V1, y = V2, color = Zoo)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = Zoo), linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  labs(
    title = "PCoA (Bray-Curtis) on pathways",
    x = paste0("Axis 1 (", eig_pct_bray[1], "%)"),
    y = paste0("Axis 2 (", eig_pct_bray[2], "%)")
  ) +
  theme_bw()

ggplot2::ggsave(
  filename = file.path(home, "results", "pcoa_pathways_bray.png"),
  plot = pcoa_bray_plot,
  width = 8,
  height = 6,
  dpi = 300
)


