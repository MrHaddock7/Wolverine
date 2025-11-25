library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

ps


filename <- file.choose()
ps_clean <- readRDS(filename)
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps_clean, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Zoo", title="Bray NMDS")
p <- plot_ordination(ps.prop, ord.nmds.bray, color="Zoo", title="Bray NMDS")

# Extract the factor levels (non-NA Zoo groups)
zlv <- levels(sample_data(ps.prop)$Zoo)

# Create a color palette for the real Zoo groups
pal <- setNames(
  RColorBrewer::brewer.pal(length(zlv), "Set2"),
  zlv
)

# Add NA = black and plot
p + scale_color_manual(
  values = pal,
  na.value = "black"
)

top40 <- names(sort(taxa_sums(ps_clean), decreasing=TRUE))[1:40]
ps.top40 <- transform_sample_counts(ps_clean, function(OTU) OTU/sum(OTU))
ps.top40 <- prune_taxa(top40, ps.top40)
plot_bar(ps.top40, x="Age", fill="Phylum") + facet_wrap(~Zoo, scales="free_x")  +
  geom_bar(stat="identity", color=NA) 

