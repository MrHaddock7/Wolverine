library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names


asv_per_sample <- apply(otu_table(ps_no_ctrl), 1, function(x) sum(x > 0))
asv_per_sample


###Check sample depth
sample_depths <- sample_sums(ps_1126)
sample_depths
bad_depth <- names(sample_depths[sample_depths < 15000])
bad_depth

###Check alpha diversity and observed ASV count per sample
library(phyloseq)
alpha <- estimate_richness(ps_no_ctrl, measures=c("Shannon", "Observed"))
alpha
write.csv(alpha, "alpha_diversity.csv", row.names = TRUE)

###Plot alpha
library(phyloseq)
library(ggplot2)
library(dplyr)

# Example: calculate Shannon diversity
alpha_div <- estimate_richness(ps_no_ctrl, measures = "Shannon")

# Combine with sample metadata
meta <- data.frame(sample_data(ps_no_ctrl))
alpha_meta <- cbind(alpha_div, meta)
alpha_meta$Sex
# Ensure Sex is a factor and include NA as a level
alpha_meta$Sex <- factor(ifelse(is.na(alpha_meta$Sex), "NA",
                                ifelse(alpha_meta$Sex == "1", "Female", 
                                       ifelse(alpha_meta$Sex == 0, "Male", "NA"))),
                         levels = c("Female", "Male", "NA"))

# Plot
ggplot(alpha_meta, aes(x = Zoo, y = Shannon, fill = Zoo)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Sex", y = "Shannon Diversity") #+
  #scale_fill_manual(values = c("Female" = "pink", "Male" = "lightblue", "NA" = "grey"))


ggsave("alpha_box_per_zoo.png", width = 6, height = 5, dpi = 300)

### Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps_1126, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

p <- plot_ordination(ps.prop, ord.nmds.bray, color="Zoo", title="Bray NMDS")

p + geom_label(aes(label = NGI.ID), 
               size = 3, 
               alpha = 0.7)

ggsave("NMDS_zoo_individual.png", width = 6, height = 5, dpi = 300)

####


top40 <- names(sort(taxa_sums(ps_no_ctrl), decreasing=TRUE))[1:40]
ps.top40 <- transform_sample_counts(ps_no_ctrl, function(OTU) OTU/sum(OTU))
ps.top40 <- prune_taxa(top40, ps.top40)
plot_bar(ps.top40, x="Age", fill="Phylum") + facet_wrap(~Zoo, scales="free_x")  +
  geom_bar(stat="identity", color=NA) 


##Plot beta diversity
library(phyloseq)
library(vegan)
library(ggplot2)   
sample_data(ps_1126)$Zoo <- as.factor(sample_data(ps_1126)$Zoo)
bc_dist <- phyloseq::distance(ps_1126, method = "bray")
ordination <- ordinate(ps_1126, method = "PCoA", distance = bc_dist)
p <- plot_ordination(ps_1126, ordination, color = "Zoo") +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Beta Diversity (PCoA) of Wolverines in Zoos")

print(p)


metadata <- as(sample_data(ps_1126), "data.frame")
adonis_result <- vegan::adonis2(bc_dist ~ Zoo, data = metadata)
adonis_result


###Adonis2
library(phyloseq)
library(vegan)

# Convert phyloseq sample data to a data.frame
meta <- data.frame(sample_data(ps_merged_1201))

# Convert relevant variables to factors
factor_vars <- c("Zoo", "Sex", "Sweden", "Meat", "Fish", "Fruit", "Individual")
meta[factor_vars] <- lapply(meta[factor_vars], factor)

# Check missing values
colSums(is.na(meta))

adonis_results <- list()
# Keep only samples with non-missing Sex
for (col in colnames(meta)) {
  if (length(unique(meta[[col]])) <= 1) next
  
  # Keep only samples without NAs
  meta_clean <- meta[complete.cases(meta[[col]]), ]
  
  if (length(unique(meta_clean[[col]])) < 2) next
  
  # Subset distance matrix
  common_samples <- intersect(rownames(meta_clean), rownames(as.matrix(bc_dist)))
  bc_dist_clean <- as.dist(as.matrix(bc_dist)[common_samples, common_samples])
  meta_clean <- meta_clean[common_samples, ]
  
  # Construct formula dynamically
  fmla <- as.formula(paste("bc_dist_clean ~", col))
  
  # Run PERMANOVA
  adonis_res <- adonis2(fmla, data = meta_clean, permutations = 999)
  
  # Store results
  adonis_results[[col]] <- adonis_res
}
adonis_results

save(adonis_results, file = "C:/Users/Lovisa/Documents/Wolverine/data/adonis_result.txt")
