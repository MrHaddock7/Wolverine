# Load libraries
install.packages("betapart")
library(phyloseq)
library(betapart)
library(ggplot2)
library(dplyr)

# Extract OTU table and metadata
otu_table_ps <- as.data.frame(otu_table(ps_rare))
meta <- as.data.frame(sample_data(ps_rare))

# Make sure OTU table has samples as rows
if(taxa_are_rows(ps_rare)){
  otu_table_ps <- t(otu_table_ps)
}

# Convert to presence/absence
otu_pa <- ifelse(otu_table_ps > 0, 1, 0)

# Compute beta diversity (Sørensen)
beta_sor <- beta.pair(otu_pa, index.family = "sorensen") 
# beta_sor$beta.sor = total, beta_sor$beta.sim = turnover, beta_sor$beta.sne = nestedness

#  Create pairwise data frame
beta_matrix <- as.matrix(beta_sor$beta.sor)
samples <- rownames(beta_matrix)

beta_df <- expand.grid(sample1 = samples,
                       sample2 = samples) %>%
  filter(sample1 != sample2)  # remove self-comparisons

# Add beta values
beta_df$beta_total <- mapply(function(x, y) beta_matrix[x, y],
                             beta_df$sample1, beta_df$sample2)

# Also add turnover and nestedness
beta_df$beta_turnover <- mapply(function(x, y) as.matrix(beta_sor$beta.sim)[x, y],
                                beta_df$sample1, beta_df$sample2)
beta_df$beta_nestedness <- mapply(function(x, y) as.matrix(beta_sor$beta.sne)[x, y],
                                  beta_df$sample1, beta_df$sample2)

# Add metadata
beta_df$Individual1 <- meta$Individual[match(beta_df$sample1, rownames(meta))]
beta_df$Individual2 <- meta$Individual[match(beta_df$sample2, rownames(meta))]
beta_df$Zoo1 <- meta$Zoo[match(beta_df$sample1, rownames(meta))]
beta_df$Zoo2 <- meta$Zoo[match(beta_df$sample2, rownames(meta))]

#  Classify comparisons
beta_df <- beta_df %>%
  mutate(within_ind = ifelse(Individual1 == Individual2, "within", "between"),
         within_zoo = ifelse(Zoo1 == Zoo2, "within_zoo", "between_zoos"))

# Statistical tests
wilcox_total_ind <- wilcox.test(beta_total ~ within_ind, data = beta_df)
wilcox_turnover_ind <- wilcox.test(beta_turnover ~ within_ind, data = beta_df)
wilcox_nested_ind <- wilcox.test(beta_nestedness ~ within_ind, data = beta_df)

wilcox_total_zoo <- wilcox.test(beta_total ~ within_zoo, data = beta_df)
wilcox_turnover_zoo <- wilcox.test(beta_turnover ~ within_zoo, data = beta_df)
wilcox_nested_zoo <- wilcox.test(beta_nestedness ~ within_zoo, data = beta_df)

# Print results
wilcox_total_ind
wilcox_turnover_ind
wilcox_nested_ind

wilcox_total_zoo
wilcox_turnover_zoo
wilcox_nested_zoo

# 8️⃣ Visualization
plot_beta <- function(df, beta_col, title){
  ggplot(df, aes(x = within_ind, y = .data[[beta_col]])) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ylab(beta_col) +
    xlab("Comparison within vs between individual") +
    ggtitle(title) +
    theme_minimal()
}

plot_beta(beta_df, "beta_total", "Beta diversity (Sørensen) - Total")
plot_beta(beta_df, "beta_turnover", "Beta diversity - Turnover")
plot_beta(beta_df, "beta_nestedness", "Beta diversity - Nestedness")

# Similarly for zoo comparison
plot_beta_zoo <- function(df, beta_col, title){
  ggplot(df, aes(x = within_zoo, y = .data[[beta_col]])) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ylab(beta_col) +
    xlab("Comparison within vs between zoos") +
    ggtitle(title) +
    theme_minimal()
}

plot_beta_zoo(beta_df, "beta_total", "Beta diversity (Sørensen) - Total (Zoo)")
plot_beta_zoo(beta_df, "beta_turnover", "Beta diversity - Turnover (Zoo)")
plot_beta_zoo(beta_df, "beta_nestedness", "Beta diversity - Nestedness (Zoo)")



# ---- NMDS PLOT OF BETAPART SØRENSEN ----

library(vegan)

# Convert Sørensen matrix to dist
dist_sor <- as.dist(beta_sor$beta.sor)

# Run NMDS
set.seed(123)
nmds_sor <- metaMDS(dist_sor, k = 2)

# Prepare dataframe
nmds_scores <- as.data.frame(scores(nmds_sor))
nmds_scores$SampleID <- rownames(nmds_scores)
nmds_scores$Individual <- meta$Individual[match(nmds_scores$SampleID, rownames(meta))]

# Plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Individual)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  ggtitle("NMDS of Sørensen Beta Diversity (betapart)")
