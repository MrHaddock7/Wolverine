# Libraries / packages

library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)  # for stat_compare_means / stat_kruskal_test

# Loading data
# Make sure to have run rarefaction first or to have access to the ps_rel file
# These paths should be pointing towards ps_rel
source("config.R")

path_ps  <- file.path(home, "data", "ps_rel_1209.RData")
save_path <- file.path(home, "plots_and_results")

load(path_ps)
ps <- ps_rel

################################################################################
######## Relative abundance box-plots ###########################################
################################################################################

# This is for the bar-plots showing relative abundance of all samples #
# I'll give example on phylum level, but it's then easy to change the taxonomic level 
# and re-do it on family or genus level for example

# Normalize
ps_norm <- microbiome::transform(ps_rel, "compositional")
barplot_phylum <- psmelt(ps_norm)

# Replacing Sample_Type with NGI.ID
# Turn sample_ID into factor, since they are numbers but don't want them to behave as numbers
# or whatever you called that column
# also turn whatever variable you're interested in comparing (for example location)
# into a factor! For me it was zoo
barplot_phylum$NGI.ID <- as.factor(barplot_phylum$NGI.ID)  
barplot_phylum$Zoo <- as.factor(barplot_phylum$Zoo)

# Top phyla
taxa_summary <- aggregate(Abundance ~ Phylum, data = barplot_phylum, sum)
taxa_summary <- taxa_summary[order(-taxa_summary$Abundance), ]
top_phyla <- as.character(taxa_summary$Phylum[1:5])
top_phyla

# For-loop and if then function that checks to see if any of the names in the list "top" is in the 
# barplot_phylum$Phylum column, leaves unchanged if it is, changes it to "Other" if it isn't
# Re-name all uncommon phyla as "Other"
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in seq_along(barplot_phylum$Phylum)) {
  if (barplot_phylum$Phylum[i] %!in% top_phyla) {
    barplot_phylum$Phylum[i] <- "Other"
  }
}

# Re-name NA:s as "Other"
barplot_phylum$Phylum[is.na(barplot_phylum$Phylum)] <- "Unassigned"
barplot_phylum$Phylum <- factor(barplot_phylum$Phylum, levels= c(top_phyla, "Other"))

# Change order so "Other" is last, change the rest to whichever order you prefer
barplot_phylum$Phylum <- factor(barplot_phylum$Phylum, levels= c("Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota","Chloroflexota", "Other"))
levels(barplot_phylum$Phylum) # check

# Colours 
top_colors <- c(
  "Bacillota"       = "#009E73",
  "Pseudomonadota"  = "#56B4E9",
  "Actinomycetota"  = "#F0E442",
  "Bacteroidota"    = "#0072B2",
  "Chloroflexota"   = "#D55E00",
  "Other"           = "#E69F00"
)


# Plot per individual sample sorted by sample type (or other variable)
# In my case I have "Sample_ID", which is the number for all my samples, 
# but your data set might have this organized differently
relab_plot <- ggplot(barplot_phylum, aes(x = NGI.ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") +
  facet_grid(~Zoo, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = top_colors)

relab_plot
# I am also using something called "facet_grid", that separates out my variable of
# interest (sample type) into different grids.

# Saving plot, remember to change date
ggsave(
  filename = file.path(save_path, "relative_abundance_plot_1209.png"),
  plot = relab_plot,
  width = 12,
  height = 6,
  units = "in"
)

########################## Merged barplot #######################################

# This might be of interest for you, if you want to merge the samples by 
# one variable, for example sex or location or in my case sample type

colnames(sample_data(ps_norm))

# Merging all samples with the same Zoo value 

sample_data(ps_norm)$Zoo <- as.factor(sample_data(ps_norm)$Zoo)
table(sample_data(ps_norm)$Zoo, useNA="ifany")

# Merge raw counts by Zoo
ps_merged_phylum <- merge_samples(ps, "Zoo")  

# Normalize after merging
merged_phylum_rel <- transform_sample_counts(ps_merged_phylum, function(x) x / sum(x))

# Melt for plotting
barplot_phylum_merged <- psmelt(merged_phylum_rel)

# Merge rare phyla into "Other"
barplot_phylum_merged$Phylum[barplot_phylum_merged$Phylum==""] <- NA
'%!in%' <- function(x,y)!('%in%'(x,y))
for(i in seq_along(barplot_phylum_merged$Phylum)) {
  if(barplot_phylum_merged$Phylum[i] %!in% top_phyla) {
    barplot_phylum_merged$Phylum[i] <- "Other"
  }
}
barplot_phylum_merged$Phylum[is.na(barplot_phylum_merged$Phylum)] <- "Other"
barplot_phylum_merged$Phylum <- factor(barplot_phylum_merged$Phylum, 
                                       levels=c("Bacillota","Pseudomonadota","Actinomycetota",
                                                "Bacteroidota","Chloroflexota","Other"))

# Plot using the Sample column already in barplot_phylum_merged
merged_relab_plot_zoo <- ggplot(barplot_phylum_merged, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  labs(x = "Zoo", y = "Relative Abundance") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = top_colors)

merged_relab_plot_zoo

# Saving plot, change dates
ggsave(
  filename = file.path(save_path, "relative_abundance_merged_by_zoo_1209.png"),
  plot = relab_plot,
  width = 12,
  height = 6,
  units = "in"
)
