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

# Laddar metadata och picrust2 abundance data
source("config.R")

metadata_path <- file.path(home, "data", "metadata.csv")
pathway_abun_path <- file.path(
  home,
  "results/picrust2_output_26_11_2025/picrust2_out_pipeline_run/pathways_out/path_abun_unstrat.tsv.gz"
)

metadata <- read.csv(metadata_path, sep = ";", fileEncoding = "latin1")

rownames(metadata) <- metadata$NGI_ID

metacyc_abundance <- readr::read_tsv(pathway_abun_path, show_col_types = FALSE)
pathway_abundance <- metacyc_abundance %>% tibble::column_to_rownames("pathway")


# Låg statistisk kraft för vissa pathways som finns i färre än 3 samples

keep <- rowSums(pathway_abundance > 0) >=3
pathway_abundance <- pathway_abundance[keep, , drop = FALSE]


# Ta bort samples av dålig kvalite
pathway_abundance <- subset(pathway_abundance, select = -c(P34104_331,P34104_333))


# Få bort orelevanta prover från metadatan

common <- intersect(colnames(pathway_abundance), rownames(metadata))
if (length(common) == 0) stop("Inga gemensamma prov mellan metadata och abundans.")

pathway_abundance <- pathway_abundance[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]


# Differential Abundance Analysis

daa_res <- pathway_daa(
  abundance  = pathway_abundance,
  metadata   = metadata,
  group      = "Zoo",
  daa_method = "LinDA"
)

daa_res_sig <- filter(daa_res, p_adjust<0.05, abs(log2FoldChange)>=1)
nrow(daa_res)
nrow(daa_res_sig)


# Pathway daa_analysis

annotated_daa_res <- pathway_daa(
  abundance = pathway_abundance,
  metadata = metadata,
  group = 'Zoo',
  daa_method = "ALDEx2"
)

head(annotated_daa_res)

annotations <- pathway_annotation(
  daa_results_df = annotated_daa_res,
  pathway = "MetaCyc",
)

head(annotation)
