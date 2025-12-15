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

tax_df <- data.frame(
  Kingdom = rownames(pathway_abundance),
  row.names = rownames((pathway_abundance)),
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
