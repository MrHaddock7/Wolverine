library(phyloseq)

# Calculate alpha diversity (Shannon index in this example)
alpha_div <- estimate_richness(ps_rare, measures = "Shannon")

# Extract metadata for individual and zoo
meta <- as.data.frame(sample_data(ps_rare))
meta_subset <- meta[, c("Individual", "Zoo")]

# Combine metadata with alpha diversity
alpha_df <- cbind(meta_subset, alpha_div)

# Rename the alpha diversity column for clarity
colnames(alpha_df)[ncol(alpha_df)] <- "alpha_diversity"

# View the result
head(alpha_df)


library(phyloseq)
library(vegan)
library(dplyr)

# Extract OTU table and metadata
otu_mat <- as(otu_table(ps_rare), "matrix")
if(taxa_are_rows(ps_rare)) otu_mat <- t(otu_mat)  # samples in rows

meta_df <- as.data.frame(sample_data(ps_rare))
meta_df$SampleID <- rownames(meta_df)

# Compute Bray-Curtis distance
bc_dist <- vegdist(otu_mat, method = "bray")

# Convert distance matrix to a matrix
bc_mat <- as.matrix(bc_dist)

# Compute mean beta diversity **within each zoo**
beta_per_zoo <- meta_df %>%
  group_by(Zoo) %>%
  summarise(beta_diversity = {
    samples <- SampleID
    sub_dist <- bc_mat[samples, samples]  # subset distance matrix for this zoo
    mean(sub_dist[upper.tri(sub_dist)])    # mean of upper triangle
  })
beta_per_zoo

# Merge back so each individual has its zoo's beta diversity
final_df <- alpha_df %>%
  left_join(beta_per_zoo, by = "Zoo")
final_df
save(final_df, file = "C:/Users/Lovisa/Documents/Wolverine/scripts/report_generation/div_df.RData")
