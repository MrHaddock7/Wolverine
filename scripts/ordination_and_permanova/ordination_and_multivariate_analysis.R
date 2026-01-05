###This script was used for the ordination and multivariate analyses

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names



### Ordination
ps.prop <- transform_sample_counts(ps_merged_1201, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

p <- plot_ordination(ps.prop, ord.nmds.bray, color="Zoo", title="Bray NMDS")

p + geom_label(aes(label = Individual), 
               size = 3, 
               alpha = 0.7)  + ggtitle(
  paste0(
    "Bray–Curtis NMDS (stress = ",
    round(ord.nmds.bray$stress, 3),
    ")"
  )
)

ggsave("NMDS_zoo_individual.png", width = 6, height = 5, dpi = 300)


###Adonis2
library(phyloseq)
library(vegan)
sample_data(ps_filtered)
# Convert phyloseq sample data to a data.frame
ps_filtered2 <- subset_samples(
  ps_filtered,
  Zoo %in% c("Vildriket", "Helsinki", "Ähtäri", "Borås")
)
sample_data(ps_filtered2)

meta <- data.frame(sample_data(ps_filtered2))
ps_rel <- phyloseq::transform_sample_counts(ps_filtered2,  function(x) x / sum(x))
bc_dist <- phyloseq::distance(ps_rel, method = "bray")

# Convert relevant variables to factors
factor_vars <- c("Zoo", "Sex", "Sweden", "Meat", "Fish", "Fruit", "Individual")
meta[factor_vars] <- lapply(meta[factor_vars], factor)
meta$Latitude <- as.numeric(gsub(",", ".", meta$Latitude))

numeric_vars <- c("Latitude", "Age", "Enclosure_size")
meta[numeric_vars] <- lapply(meta[numeric_vars], function(x) as.numeric(as.character(x)))

# Check missing values
colSums(is.na(meta))

adonis_results <- list()
# Keep only samples with non-missing Sex
for (col in colnames(meta)) {
  # Remove NAs for this variable
  meta_clean <- meta[!is.na(meta[[col]]), ]
  
  # Skip if less than 2 unique values remain
  if (length(unique(meta_clean[[col]])) < 2) next
  
  # Subset distance matrix
  common_samples <- intersect(rownames(meta_clean), rownames(as.matrix(bc_dist)))
  bc_dist_clean <- as.dist(as.matrix(bc_dist)[common_samples, common_samples])
  meta_clean <- meta_clean[common_samples, ]
  
  # Skip if less than 2 samples remain
  if (nrow(meta_clean) < 2) next
  
  # Construct formula
  fmla <- as.formula(paste("bc_dist_clean ~", col))
  
  # Run PERMANOVA
  adonis_res <- adonis2(fmla, data = meta_clean, permutations = 999)
  
  # Store result
  adonis_results[[col]] <- adonis_res
}
adonis_results
save(adonis_results, file = "C:/Users/Lovisa/Documents/Wolverine/data/adonis_result.RData")

###Multivariate ADONIS
# Define variables in the model
model_vars <- c("Zoo", "Age", "Sex", "Latitude")

# Keep only rows with no NAs in these variables
meta_clean2 <- meta_clean[complete.cases(meta_clean[, model_vars]), ]

# Subset distance matrix to match these samples
common_samples <- intersect(rownames(meta_clean2), rownames(as.matrix(bc_dist)))
bc_dist_clean2 <- as.dist(as.matrix(bc_dist)[common_samples, common_samples])
meta_clean2 <- meta_clean2[common_samples, ]

# Run multivariable PERMANOVA
adonis_res2 <- adonis2(
  bc_dist_clean2 ~ Sex + Zoo + Age + Latitude,
  data = meta_clean2,
  permutations = 999,
  by = "margin"
)

adonis_res2



###Pairwise adonis
install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(phyloseq)
library(vegan)
library(pairwiseAdonis)

# run pairwise adonis

ps_rel <- phyloseq::transform_sample_counts(ps_filtered,  function(x) x / sum(x))
bc_dist <- phyloseq::distance(ps_rel, method = "bray")

pw <- pairwise.adonis(bc_dist, phyloseq::sample_data(ps_rel)$Zoo)
print(pw)





