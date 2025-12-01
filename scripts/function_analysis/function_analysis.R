library(ggpicrust2)
library(phyloseq)
library(dplyr)
library(readr)
library(MicrobiomeStat)
library(ggplot2)
library(tibble)
source("config.R")


metadata <- read.csv(paste(home, 'data/metadata.csv', sep = ''), sep = ';', fileEncoding = "latin1")
rownames(metadata) <- metadata$NGI.ID

metacyc_abundance <- read_tsv(
  paste(home, "results/picrust2_output_26_11_2025/picrust2_out_pipeline_run/pathways_out/path_abun_unstrat.tsv.gz", sep = '')
)

# Första kolumnen brukar heta "pathway" – gör den till rownames:
pathway_abundance <- metacyc_abundance %>%
  tibble::column_to_rownames("pathway")

head(metadata)
head(pathway_abundance)


## 3. Synka metadata och abundans (bara gemensamma prover)
common_samples <- intersect(colnames(pathway_abundance), rownames(metadata))

length(common_samples)


pathway_abundance <- pathway_abundance[, common_samples]
metadata <- metadata[common_samples, ]

## 4. Välj en gruppvariabel – t.ex. Zoo
table(metadata$Zoo)

# Om det ser ut som Skansen / Helsinki / Gaia / ... med flera prover per grupp
# kan vi köra DAA så här:

daa_res <- pathway_daa(
  abundance  = pathway_abundance,
  metadata   = metadata,
  group      = "Zoo",   # <-- VIKTIGT: en riktig grupp-kolumn, inte NGI.ID
  daa_method = "LinDA"
)

## 5. Visualisera resultat
visualize_daa(daa_res)








keep <- rowSums(pathway_abundance > 0) >= 3  # eller >= 5 om du vill vara strängare
sum(keep)  # hur många pathways överlever?

pathway_abundance_filt <- pathway_abundance[keep, ]

head(pathway_abundance_filt)

daa_res <- pathway_daa(
  abundance  = pathway_abundance_filt,
  metadata   = metadata,
  group      = "Zoo",
  daa_method = "LinDA"
)


head(daa_res$daa_results_df)
colnames(daa_res$daa_results_df)





path_tbl <- pathway_abundance_filt %>%
  rownames_to_column("pathway")

top20 <- path_tbl %>%
  mutate(total = rowSums(across(-pathway))) %>%
  arrange(desc(total)) %>%
  slice_head(n = 20)

pathway_abundance_top20 <- top20 %>%
  select(-total) %>%
  column_to_rownames("pathway")



pathway_heatmap(
  abundance = pathway_abundance_top20,
  metadata = metadata,
  group = "Zoo",

) +
  labs(title = "Pathway Abundance Heatmap")




