library(phyloseq)
#Merge ps 
sample_names(ps_rare)
sample_data(ps_rare)["P34104_257", "Individual"] <- "SK-U"

names_df <- read.csv("C:/Users/Lovisa/Documents/Wolverine/data/NGI_ID_to_Name.csv", sep = ";", stringsAsFactors = FALSE)


ps_filt <- prune_samples(!(sample_names(ps_rare) %in% c("P34104_331", "P34104_333")), ps_rare)
sample_names(ps_filt)
ps_merged <- merge_samples(ps_filt, "Individual")
meta <- data.frame(sample_data(ps_filt))
meta$Name <- names_df$Name[match(rownames(meta), names_df$NGI.ID)]
meta_merged <- meta[!duplicated(meta$Individual), ]
rownames(meta_merged) <- meta_merged$Individual
sample_data(ps_merged) <- sample_data(meta_merged)

# Calculate alpha diversity (Shannon index in this example)
alpha_div <- estimate_richness(ps_merged, measures = "Shannon")

# Extract metadata for individual and zoo
meta <- as.data.frame(sample_data(ps_merged))
meta_subset <- meta[, c("Individual", "Zoo", "Name")]

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
