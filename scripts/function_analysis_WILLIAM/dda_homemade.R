library(microeco)
library(file2meco)
library(dplyr)
library(tibble)
library(readr)
library(ggpicrust2)
library(phyloseq)
library(ggpubr)

source("config.R")

metadata <- read.csv(file.path(home, "data", "metadata.csv"), sep = ";", fileEncoding = "UTF-8")
rownames(metadata) <- metadata$NGI.ID

pathway_out <- readr::read_tsv("results/picrust2_output_26_11_2025/picrust2_out_pipeline_run/pathways_out/path_abun_unstrat.tsv")
ref <- ggpicrust2:::load_reference_data("MetaCyc")


tax_df <- data.frame(
  Kingdom = df1$description,
  row.names = df1$pathway,
  stringsAsFactors = FALSE
)

tax <- tax_table(as.matrix(tax_df)) 

samples <- otu_table(as.matrix(pathway_abundance), taxa_are_rows = TRUE)

ps_pathways <- phyloseq(samples, tax_table = tax, sam_data = sample_data(metadata))



data <- phyloseq2meco(ps_pathways)

#calculate relative abundances
data$cal_abund()

#now calculate the biodindecator taxa with differential abundance using random forest (rf) at bacterial genus level
t2 <- trans_diff$new(dataset = data, method = "KW", group = "Zoo", taxa_level = "Kingdom")


t2$plot_diff_abund(order_x_mean = TRUE) ##plot

res_wilcox <- t2$res_diff ##table
