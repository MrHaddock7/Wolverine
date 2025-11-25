# Install DADA2 package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2") #change version if updated

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")

# Follow instructions to install cutadapt, you will do that on your computer and then
# run it remotely in R, https://cutadapt.readthedocs.io/en/stable/installation.html

##### 

# load packages, check which version you have 

library(dada2); packageVersion("dada2")
library(ShortRead)
library(cutadapt)
library(ggplot2)
library(MASS)

library(phyloseq)
library(Biostrings)
library(decontam)
library(microbiome)
library(microbiomeutilities)
library(vegan)
library(eulerr)

# Some extra packages to use for plotting, should be easy to find and directly install
library(ggplot2)
library(tidyverse)
library(grid)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(ggsignif)

# you may need to change what kind of data separation you have in your metadata file
### 1. Load metadata ----------------------------------------------------------

metadata <- read.csv(
  "C:/Users/anton/Desktop/Applied Bioinformatics/metadata.csv",
  sep = ";",
  stringsAsFactors = FALSE
)

# Set rownames to NGI.ID
rownames(metadata) <- metadata$NGI.ID


### 2. Extract sample names from seqtab.nochim --------------------------------

samples.out <- rownames(seqtab.nochim)


### 3. Clean sample names to match metadata -----------------------------------

# Remove trailing underscore
samples.clean <- sub("_$", "", samples.out)

# Remove Illumina sample index (_S###)
samples.clean <- sub("_S[0-9]+$", "", samples.clean)

# Now assign cleaned names to seqtab.nochim
rownames(seqtab.nochim) <- samples.clean


### 4. Filter metadata to only samples that have ASV data ----------------------

metadata.filtered <- metadata[samples.clean, ]

# Ensure rownames match
rownames(metadata.filtered) <- samples.clean


### 5. Build phyloseq object --------------------------------------------------

# Metadata must NOT be transposed
t_metadata <- metadata.filtered

ps_raw <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(t_metadata),
  tax_table(taxa)
)


### 6. Convert sequences into ASV names ---------------------------------------

dna <- Biostrings::DNAStringSet(taxa_names(ps_raw))
names(dna) <- taxa_names(ps_raw)
ps_raw <- merge_phyloseq(ps_raw, dna)

taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))


### 7. Clean up workspace (optional) ------------------------------------------

rm(
  dadaFs, dadaRs, errF, errR, mergers, out_trim, seqtab, seqtab2, track2,
  filtFs, filtRs, fnFs, fnRs, FWD, FWD.orients, path, REV, REV.orients
)

ps_raw

####Cleaning PS#####
get_taxa_unique(ps_raw, "Family")
get_taxa_unique(ps_raw, "Class")
get_taxa_unique(ps_raw, "Order")

mitochondria <- phyloseq::subset_taxa(ps_raw, Family=="Mitochondria")
mitochondria
sample_sums(mitochondria)[order(sample_sums(mitochondria))]

chloroplast <- subset_taxa(ps_raw, Order=="Chloroplast")
chloroplast
tax_table(chloroplast)[,1:7]
sample_sums(chloroplast)[order(sample_sums(chloroplast))]

colnames(tax_table(ps_raw))
table(tax_table(ps_raw)[, "Kingdom"], useNA = "always")


# Eukaryotes and unassigned (I don't have Eukaryotes in my samples)
unassigned <- subset_taxa(ps_raw, is.na(Kingdom))
unassigned
tax_table(unassigned)[,1:7]
sample_sums(unassigned)[order(sample_sums(unassigned))]

#Remove eukaryotes
eukaryotes <- subset_taxa(ps_raw, Kingdom == "Eukaryota")

#Remove archaea
archaea <- subset_taxa(ps_raw, Kingdom == "Archaea")
tax_table(eukaryotes)[, 1:7]
tax_table(archaea)[, 1:7]
sample_sums(eukaryotes)[order(sample_sums(eukaryotes))]
sample_sums(archaea)[order(sample_sums(archaea))]


# filter mitochondria 
badTaxa <- row.names(tax_table(mitochondria)[,1])
allTaxa <- taxa_names(ps_raw)
goodTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_raw <- prune_taxa(goodTaxa, ps_raw)
# filter chloroplast
badTaxa <- row.names(tax_table(chloroplast)[,1])
allTaxa <- taxa_names(ps_raw)
goodTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_raw <- prune_taxa(goodTaxa, ps_raw)
# filter NA reads
badTaxa <- row.names(tax_table(unassigned)[,1])
allTaxa <- taxa_names(ps_raw)
goodTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_raw <- prune_taxa(goodTaxa, ps_raw)
# filter eukaryotes reads
badTaxa <- row.names(tax_table(eukaryotes)[,1])
allTaxa <- taxa_names(ps_raw)
goodTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_raw <- prune_taxa(goodTaxa, ps_raw)
# filter archaea reads
badTaxa <- row.names(tax_table(archaea)[,1])
allTaxa <- taxa_names(ps_raw)
goodTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_raw <- prune_taxa(goodTaxa, ps_raw)


# filter zero ASVs
table(taxa_sums(ps_raw) < 1)
ps_raw
ps_raw <- prune_taxa(taxa_sums(ps_raw) > 0, ps_raw)
ps_raw
# clean up some objects in work-space
rm(mitochondria, chloroplast, allTaxa, goodTaxa, badTaxa)

################################################################################
########################## Decontam ############################################
################################################################################

# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# the website has a very nice introduction that we will follow

# We will likely not have the DNA concentration for each sample unless you performed
# all lab work yourself, in which case you can add the frequency method to your
# script. It is however not necessary and we will use the prevalence method

# The prevalence method requires a column in your metadata that says whether your
# sample is a true sample or control. 

# Check metadata
head(sample_data(ps_raw))

# Convert sample_data to data.frame for plotting
df <- as.data.frame(sample_data(ps_raw))
df$LibrarySize <- sample_sums(ps_raw)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#v Converting sample data to data frame with correct ASCII
# Ensure User_ID is character
df$User_ID <- as.character(df$User_ID)

# Convert to UTF-8 to avoid translation warnings
df$User_ID <- iconv(df$User_ID, from = "", to = "UTF-8")

# Now classify control vs sample
df$Sample_or_Control <- ifelse(grepl("^Control_", df$User_ID),
                               "Control",
                               "Sample")


# Create Sample_or_Control column based on User_ID
df$Sample_or_Control <- ifelse(grepl("^Control_", df$User_ID),
                               "Control",
                               "Sample")

# Plot library sizes colored by Sample_or_Control
ggplot(df, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) +
  geom_point()

# Add Sample_or_Control to ps_raw sample_data
sample_data(ps_raw)$Sample_or_Control <- df$Sample_or_Control[match(
  rownames(sample_data(ps_raw)), df$NGI.ID
)]

# Add negative control indicator
sample_data(ps_raw)$is.neg <- sample_data(ps_raw)$Sample_or_Control == "Control"

# Identify contaminants using prevalence method (threshold = 0.2)
contamdf.prev <- isContaminant(ps_raw, method="prevalence", neg="is.neg", threshold=0.2)
table(contamdf.prev$contaminant)



# Transform to presence/absence for plotting
ps.pa <- transform_sample_counts(ps_raw, function(abund) 1*(abund > 0))

# Test
# Transform counts to presence/absence for prevalence
ps.pa <- transform_sample_counts(ps_raw, function(x) 1*(x>0))

# Split into negative controls and true samples
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)

# Calculate prevalence for each ASV
prev.df <- data.frame(
  ASV = taxa_names(ps.pa),
  pa.neg = taxa_sums(ps.pa.neg) / nsamples(ps.pa.neg),  # prevalence in controls
  pa.pos = taxa_sums(ps.pa.pos) / nsamples(ps.pa.pos)   # prevalence in true samples
)

# Run decontam with your chosen threshold
threshold <- 0.2  # adjust as needed
contamdf.prev <- isContaminant(ps_raw, method="prevalence", neg="is.neg", threshold=threshold)

# Add contaminant info
prev.df$contaminant <- contamdf.prev$contaminant

# Plot prevalence
library(ggplot2)
ggplot(prev.df, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point(size=2, alpha=0.7) +
  xlab("Prevalence in negative controls") +
  ylab("Prevalence in true samples") +
  ggtitle(paste("Contaminant detection (threshold =", threshold, ")")) +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color="gray") +
  theme(legend.title = element_blank())


# Subset to negative controls and true samples
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)

# Data frame of prevalence in positive vs negative samples
df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos),
  pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev$contaminant
)

# Plot prevalence
ggplot(data=df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)")

# Prune contaminants
ps_noncontam <- prune_taxa(!contamdf.prev$contaminant, ps_raw)

# Remove the 'is.neg' column from sample_data
sample_data(ps_noncontam)$is.neg <- NULL

# Check final sample_data
sample_data(ps_noncontam)

# Clean up workspace
rm(df, df.pa, ps.pa, ps.pa.neg, ps.pa.pos, contamdf.prev)


table(sample_data(ps_raw)$Sample_or_Control)


###############################################################################
################# More sample clean-up and preparation #########################
################################################################################

# Remove all the control samples since we don't need them anymore
# I named it after my study species, call it whatever you want.
ps_no_ctrl <- subset_samples(ps_noncontam, Sample_or_Control %in% "Sample")

##### Remove very poor samples (i.e the ones with very low number of reads)

sample_sums(ps_no_ctrl)[order(-sample_sums(ps_no_ctrl))]
ps_no_ctrl = prune_samples(sample_sums(ps_no_ctrl)>=1000, ps_no_ctrl)

# Check ps object
summarize_phyloseq(ps_no_ctrl)

n_raw  <- ntaxa(ps_raw)
n_clean <- ntaxa(ps_no_ctrl)
n_raw
n_clean
removed <- n_raw - n_clean
removed

removed / n_raw * 100

ps.pa <- transform_sample_counts(ps_raw, function(x) 1*(x>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg, ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$is.neg, ps.pa)

df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos),
  pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev$contaminant
)

ggplot(df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence in controls") +
  ylab("Prevalence in samples")

table(contamdf.prev$contaminant)

sum_before <- sum(sample_sums(ps_raw))
sum_after  <- sum(sample_sums(ps_no_ctrl))
removed_reads <- sum_before - sum_after
removed_reads
c(sum_before, sum_after, removed_reads)


# Choose a path and filename
saveRDS(ps_no_ctrl, file = "C:/Users/anton/Desktop/Applied Bioinformatics/ps_no_ctrl_clean.rds")



################################################################################
######################### Plotting and data-exploration ########################
################################################################################

load("~/path/to/Phyloseq_decontam")

# Set color-palette for plotting, you can pick your own colors here, whatever
# suits you!

color_list <- c("#8A9A5B", "#a96b17", "#d79235", "#c45c04", "#452a05", "#646444", "#353b24", "#94842b", "#302a04", "#cc6350", "#242b10")

################################################################################
######## Relative abundance box-plots ###########################################
################################################################################

#Checking before and after (relative abundance)
ps_rel_raw <- transform_sample_counts(ps_raw, function(x) x / sum(x))
ps_rel_clean <- transform_sample_counts(ps_no_ctrl, function(x) x / sum(x))

df_raw <- psmelt(ps_rel_raw)
df_clean <- psmelt(ps_rel_clean)

df_raw$Phylum[df_raw$Phylum %in% c("", " ", NA, "NA", "Unassigned")] <- "Unassigned"
df_clean$Phylum[df_clean$Phylum %in% c("", " ", NA, "NA", "Unassigned")] <- "Unassigned"

phylum_levels <- sort(unique(c(df_raw$Phylum, df_clean$Phylum)))
df_raw$Phylum <- factor(df_raw$Phylum, levels = phylum_levels)
df_clean$Phylum <- factor(df_clean$Phylum, levels = phylum_levels)

# Before pruning
ggplot(df_raw, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  ggtitle("Relative abundance BEFORE contaminant removal") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

# After pruning
ggplot(df_clean, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  ggtitle("Relative abundance AFTER contaminant removal") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

# Checking abundance
# Melt phyloseq objects
df_raw <- psmelt(ps_raw)        # before pruning
df_clean <- psmelt(ps_no_ctrl)  # after pruning

# Fix missing phylum names
df_raw$Phylum[df_raw$Phylum %in% c("", " ", NA, "NA", "Unassigned")] <- "Unassigned"
df_clean$Phylum[df_clean$Phylum %in% c("", " ", NA, "NA", "Unassigned")] <- "Unassigned"

# Plot absolute abundance before
ggplot(df_raw, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  ggtitle("Absolute abundance BEFORE pruning") +
  ylab("Read counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

# Plot absolute abundance after
ggplot(df_clean, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  ggtitle("Absolute abundance AFTER pruning") +
  ylab("Read counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))


# This is for the bar-plots showing relative abundance of all samples #
# I'll give example on phylum level, but it's then easy to change the taxonomic level 
# and re-do it on family or genus level for example

# Use the ps object that has been normalized!
barplot_phylum <- psmelt(ps_no_ctrl)
table(barplot_phylum$Phylum, useNA="ifany") # table of phyla
barplot_phylum$Phylum[barplot_phylum$Phylum==""] <- "NA" # change if there are empty names to NA
table(barplot_phylum$Phylum, useNA="ifany") # check

# calculate sum of abundances
taxa_summary <- aggregate(barplot_phylum$Abundance, by=list(Category= barplot_phylum$Phylum), FUN=sum)
taxa_summary <- taxa_summary[order(-taxa_summary$x) ,]
taxa_summary

#Print off all top 5 phyla (choose your own cut-off)
list(as.character(taxa_summary[c(1:5),1]))

top_phyla <- as.character(c("Actinobacteriota","Firmicutes","Bacteroidota","Proteobacteria","Desulfobacterota"))

# Function
'%!in%' <- function(x,y)!('%in%'(x,y))

# For-loop and if then function that checks to see if any of the names in the list "top" is in the 
# barplot_phylum$Phylum column, leaves unchanged if it is, changes it to "Other" if it isn't
# Re-name all uncommon phyla as "Other"
for (i in seq_along(barplot_phylum$Phylum)) {
  if (barplot_phylum$Phylum[i] %!in% top_phyla) {
    barplot_phylum$Phylum[i] <- "Other"
  }
}

# Re-name NA:s as "Other"
ps_no_ctrl$Phylum[is.na(barplot_phylum$Phylum)] <- "Other"

tax <- tax_table(ps_no_ctrl)
tax <- as.data.frame(tax)
tax$Phylum[
  is.na(tax$Phylum) |
    tax$Phylum == "" |
    tax$Phylum == " " |
    tax$Phylum == "NA" |
    tax$Phylum == "Unassigned"
] <- "Unassigned"
tax_table(ps_no_ctrl) <- as.matrix(tax)
barplot_phylum <- psmelt(ps_no_ctrl)
table(barplot_phylum$Phylum, useNA="ifany")
barplot_phylum


# Change order so "Other" is last, change the rest to whichever order you prefer
barplot_phylum$Phylum <- factor(barplot_phylum$Phylum, levels= c("Actinobacteriota", "Firmicutes", "Proteobacteria", "Bacteroidota", "Desulfobacterota", "Other"))
levels(barplot_phylum$Phylum) # check

# Turn sample_ID into factor, since they are numbers but don't want them to behave as numbers
# or whatever you called that column
barplot_phylum$Sample_ID <- as.factor(barplot_phylum$Sample_ID)

# also turn whatever variable you're interested in comparing (for example location)
# into a factor! For me it was sample type
barplot_phylum$Sample_type <- as.factor(barplot_phylum$Sample_type)

# Plot per individual sample sorted by sample type (or other variable)
# In my case I have "Sample_ID", which is the number for all my samples, 
# but your data set might have this organized differently

ggplot(barplot_phylum, aes(x = Sample_ID, y=Abundance, order = Phylum)) + 
  geom_bar((aes(fill=Phylum)), stat="identity", position="stack", width = 1) +
  ylab("Relative abundance") + xlab("Sample ID") + theme_bw() +
  # coord_flip() + # comment out this line if you want vertical bars
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.title.x = element_text(size=13), 
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(size=12, color ="black"), 
        axis.text.x = element_text(size=6)) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.key.size=unit(0.7,"cm") ) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=(color_list[c(10,6,8,5,4,3)])) + 
  facet_grid(~Sample_type, scales="free", space="free") + # change this per your liking!
  scale_y_continuous(expand = c(0,0)) +
  theme(strip.text.x = element_text(size = 8, color = "black"), 
        strip.background = element_rect(fill = "#ffffff"))

# This is based in ggplot2. I would recommend looking up various ways of editing
# the aesthetics of the plot, like text font, size, color etc. 
# I am also using something called "facet_grid", that separates out my variable of
# interest (sample type) into different grids.

########################## Merged barplot #######################################

# This might be of interest for you, if you want to merge the samples by 
# one variable, for example sex or location or in my case sample type

ps_merged_phylum = merge_samples(ps_Lagilis, "Sample_type")

# Rel. abundance transformation
merged_phylum_rel = transform_sample_counts(ps_merged_phylum, function(OTU) OTU/sum(OTU))

# Now merge rare phyla into "Other" as we did before
barplot_phylum_merged <- psmelt(merged_phylum_rel)
table(barplot_phylum_merged$Phylum, useNA="ifany") # table of phyla
barplot_phylum_merged$Phylum[barplot_phylum_merged$Phylum==""] <- "NA" # change if there are empty names to NA
table(barplot_phylum_merged$Phylum, useNA="ifany") # check

# Calculate sum of abundances
taxa_summary2 <- aggregate(barplot_phylum_merged$Abundance, by=list(Category = barplot_phylum_merged$Phylum), FUN=sum)
taxa_summary2 <- taxa_summary[order(-taxa_summary$x) ,]
taxa_summary2

# Choose your own cut-off threshold. Copy all names for next revalue step.
list(as.character(taxa_summary2[c(1:5),1]))

# Make list of top 5 most abundant phyla, this will probably be the same as before
top_phyla <- as.character(c("Actinobacteriota","Firmicutes","Bacteroidota","Proteobacteria","Desulfobacterota"))

# Same function
'%!in%' <- function(x,y)!('%in%'(x,y))

for (i in seq_along(barplot_phylum_merged$Phylum)) {
  if (barplot_phylum_merged$Phylum[i] %!in% top_phyla) {
    barplot_phylum_merged$Phylum[i] <- "Other"
  }
}

barplot_phylum_merged$Phylum[is.na(barplot_phylum_merged$Phylum)] <- "Other"

# change order so Other is last, change to whichever order you prefer
barplot_phylum_merged$Phylum <- factor(barplot_phylum_merged$Phylum, levels= c("Actinobacteriota", "Firmicutes", "Proteobacteria", "Bacteroidota", "Desulfobacterota", "Other"))
levels(barplot_phylum_merged$Phylum)

# Have to turn sample_ID into factor or it won't plot right
barplot_phylum_merged$Sample_ID <- as.factor(barplot_phylum_merged$Sample_ID)

# Once again I'm plotting per sample type, so you need to change it to 
# whatever other variable of interest you have
ggplot(barplot_phylum_merged, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") + theme_light() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Sample type", y = "Relative Abundance") +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(size=12, color ="black"), 
        axis.text.x = element_text(size=10)) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.key.size=unit(0.5,"cm")) +
  scale_fill_manual(values=color_list[c(10,6,8,5,4,3)])

################################################################################
###### Plot different diversity indices ########################################
################################################################################

# Make box-plots of various diversity indices, here I compare sample type but 
# change it to whatever else you want. The plots are somewhat ugly so if you
# want to use it you may need to edit it more.

Sample_type_div = plot_richness(ps_rare, x = "Sample_type", measures = c("Observed","Shannon", "Simpson"), color = "Sample_type", scales = "free_y") +
  geom_boxplot(width = 0.5) +
  scale_color_manual(values = color_list[c(2,5,8)]) +
  theme_bw() + 
  theme(legend.position = "none") +
  scale_x_discrete() +
  theme(axis.title.y = element_text(size=8), 
        axis.text.y = element_text(size=10)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())

################################################################################
########## Box-plots for ASV richness ##########################################
################################################################################

# Since we're using alpha diversity it's good to use the rarefied object
Sampletype_boxplot <- plot_richness(ps_rare, x = "Sample_type", measures = "Observed")

ggplot(Sampletype_boxplot$data, aes(x = Sample_type, y = value, fill = Sample_type)) +
  geom_boxplot(alpha = 1, color = "black", outlier.color = "white", outlier.size = 0, 
               lwd=0.36, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.25, color = "black", pch = 21) +  
  ylab("Richness") + xlab("Sample type") + theme_bw() +
  theme(axis.title.x = element_text(size=18, margin=margin(t=7)), 
        axis.title.y = element_text(size=18, margin=margin(r=8)),
        axis.text.y  = element_text(size=18, color="black"), 
        axis.text.x  = element_text(size=18, color="black"),
        strip.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "grey29"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_fill_manual(values = color_list[c(6,4,8)]) +
  scale_color_manual(values = color_list) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.key.size = unit(1,"cm")) +
  theme(legend.position = "none") +
  stat_kruskal_test(group.by = "x.var", label = "p.signif", label.y = 220) +
  # Note here that I add in a statistical test to indicate significance
  # It depends on your data which test you should use
  ggtitle("Boxplot showing difference in ASV richness between sample-types")

# Example on how you might plot two different box-plots in the same plot
# using ggarrange. You can't use par(mfrow = c(2,2)) with ggplot 
ggarrange(plot1_name, plot2_name, nrow = 1, labels = c("A", "B"))

################################################################################
################################# PCoA/PCA plot ################################
################################################################################

# Using a PCoA (principal coordinate analysis) or PCA (principal component analysis)
# is a nice way to visualize the ASV composition of your samples, i.e the beta diversity

# I ended up using Aitchinson's distance but I added examples for both
# Bray-Curtis ordination and Aitchinson's distance.
# It might be good to try different methods.

# Turn your sample variables into factors
# Important that you use the normalized data set:
sample_data(ps_rel)$Sex <- as.factor(sample_data(ps_rel)$Sex)
sample_data(ps_rel)$Sample_type <- as.factor(sample_data(ps_rel)$Sample_type)

################################################################################
# Elin's PcOA plot example

# change here the colors to what you like or want to use. Use as many colors as you have categories. 
col_gir[2] <- "#c45c1a"
col_gir[1] <- "#F2C188"
names(col_gir) <- c("Northern", "Reticulated", "Masai")

ord_ps = ordinate(psF.rel, "PCoA", "bray")
g <- plot_ordination(psF.rel, ord_ps, justDF=T)
head(g)

# extract % variation
tempx <- as.character(round(ord_ps$values$Relative_eig[1],3) *100)
tempy <- as.character(round(ord_ps$values$Relative_eig[2],3) *100)
set.seed(386)

# pcoa plot Giraffe Species
ggplot(g, aes(x=Axis.1, y=Axis.2, color=Species, fill=Species)) + theme_bw() + 
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, level=0.9)  +  # confidence interval
  geom_point(color="black", shape=21, size = 2.5, stroke=0.3) +
  xlab(paste0("PC1 (",temps,"%)")) + ylab(paste0("PC2 (",tempy,"%)")) +
  theme(axis.title.x = element_text(size=13, margin = margin(t=3)),
        axis.title.y = element_text(size=13, margin = margin(r=2)),
        axis.text.y  = element_text(size=11, color = "black"), 
        axis.text.x  = element_text(size=11, color = "black"),
        strip.text   = element_text(size=13)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 13), 
        legend.text  = element_text(size = 13),
        legend.key.size = unit(0.6, "cm"),
        legend.position = "right") +
  scale_color_manual(values = col_gir) +
  scale_fill_manual(values  = col_gir)

################################################################################
# PCA plot example

# This transforms the relative abundances using "clr" (centered log ratio) transformation
ps_clr <- microbiome::transform(ps_rel, "clr")

# Ordinate RDA (redundancy analysis), this is for the distances for the PCA
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

# calculate the relative eigenvalues for pc1 and pc2
# $CA is the object in ord_clr that includes the eigenvalues for each PC
sum1 <- sum(ord_clr$CA$eig)

tempx <- as.character(round((ord_clr$CA$eig[1]/sum1)*100))
tempy <- as.character(round((ord_clr$CA$eig[2]/sum1)*100))

PCA_aitch <- plot_ordination(ps_rel, ord_clr, justDF=T)

ggplot(PCA_aitch, aes(x=PC1, y=PC2, fill = Sample_type, color = Sample_type, group = Sample_type)) + 
  theme_bw() + 
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, level=0.95) +  # confidence interval
  geom_point(aes(shape = Sex), size = 2.5, stroke=0.3) +  # move shape aesthetic here
  labs(x=(paste0("PC1 (",tempx,"%)")), y=(paste0("PC2 (",tempy,"%)"))) +
  theme(axis.title.x = element_text(size=13, margin = margin(t=3)),
        axis.title.y = element_text(size=13, margin = margin(r=2)),
        axis.text.y  = element_text(size=11, color = "black"), 
        axis.text.x  = element_text(size=11, color = "black"),
        strip.text   = element_text(size=13)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 13), 
        legend.text  = element_text(size = 13),
        legend.key.size = unit(0.6, "cm"),
        legend.position = "right") +
  scale_color_manual(values = color_list[c(9,5,1)]) +
  scale_fill_manual(values = color_list[c(9,5,2)], guide = "none") +
  scale_shape_manual(values = c(21,22)) +  # specify shape values
  guides(color = guide_legend("Sample type")) # this is to make the legend colors work

################################################################################
############## Statistical tests on alpha diversity ############################
################################################################################

# These are tests to statistically compare the alpha diversity of different variables in the
# data. You might need to do different tests depending on your data. 

# It is probably best to use non-parametric tests with microbiome data, although you 
# can test for normality by plotting Q-Q plots, doing a quick histogram or a shapiro test. 
# A p-value > 0.05 means it's normally distributed

# hist(variable)
# qqplot(variable)
# qqline(variable)
# shapiro.test(variable)

# Get ASV richness data, you can use the "estimate_richness" function
# you could also use a different alpha diversity measure than observed richness
rich = estimate_richness(ps_rare)

# Here you can go back and use your metadata, as it will be nicely organized
# and contain all the information you'll need
# you will need to remove all the actual samples you removed before, e.g.
# your control samples and anything else you removed. 

alpha_div = subset(metadata, Sample_or_control == "Sample")

# check to see that row numbers are the same
nrow(alpha_div)
nrow(rich)

# I added an extra column with the observed richness values, you will
# need to check how many columns you have 
ncol(alpha_div)
alpha_div[,15] <- rich$Observed
names(alpha_div)[15] <- paste("Richness")

# Here I have an example of a Kruskal-Wallice test, which is the non-parametric
# one-way ANOVA, testing for two or more groups on a continuous variable

# Here I test sample types against each other, but of course change the variable
# per your data set
kw1 <- kruskal.test(Richness ~ Sample_type, data = alpha_div)

# Plot distribution, you can transform data if it's not normal, or change
# the distribution in the glm later on.
hist(alpha_div$Richness)

# And here I test a general linear model. Change the family to whatever distribution
# you have. It could for example be gaussian, poisson etc.
glm_1 <- glm(Richness ~ Sex + Temp_group, data = alpha_div, family = gaussian)
summary.glm(glm_1)

################################################################################
########## Adonis tests on ASV composition #####################################
################################################################################

# The Adonis test is included in the microbiome package, and is a version of a 
# PERMANOVA (Permutational analysis of variance). It is essentially a test for
# difference in beta diversity between sample variables. 

# Begin by using your normalized data (since we're testing beta diversity),
# and all your variables need to be as factors, so transform them as so
sample_data(ps_rel)$Sex <- as.factor(sample_data(ps_rel)$Sex)

############################## Using CLR #######################################

# centered log-ratio (CLR) transformation of counts
ps_clr <- microbiome::transform(ps_rel, "clr")

#Generate distance matrix, you can do this in various ways. 
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean")

#adonis test, on sample type, sex, location etc, whatever variable you want
# you will get a P-value, and an R2 value (how much variation is explained)
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Sex, permutations = 999)

######################## Using bray-curtis #####################################

#Generate distance matrix, you can do this in various ways. 
bray_dist_matrix <- phyloseq::distance(ps_rel, method = "bray")

#adonis test, on sample type, sex, location etc, whatever variable you want
# you will get a P-value, and an R2 value (how much variation is explained)
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps_rel)$Sex, permutations = 999)

# You can also run it on two variables, if so add a "+" between them
# if you use * you also test the interaction, which is more complicated to interpret

################################################################################
############################ Core microbiome ################################### 
################################################################################

# https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_heatmap
# https://joey711.github.io/phyloseq/plot_heatmap-examples.html

# The core microbiome can be defined in different ways. It is also somewhat
# ambiguous what the biological meaning actually is. It can however be useful
# if you want to compare your data directly to someone else's study, so I have
# included it here.

# First let us just create a regular heatmap plot of the most abundant taxa. You
# may want this to show the composition, so it's not actually the core microbiota
# yet. 

# Filter out 50 most prevalent ASVs (or however many you want)
ps_rel_pruned <- prune_taxa(names(sort(taxa_sums(ps_rel),TRUE)[1:50]), ps_rel)

tax_table(ps_rel_pruned) # get tax table for the data set

# Depending on if you have a lot of NA:s in your plot, you may want to re-name
# them to something nicer, here I have an example:

print(tax_table(ps_rel_pruned)[,5]) # this shows all top ASVs at the genus level 

# if any of the ASVs have NA:s note them down

# find the name of the order the unknown family/genus belongs to 
tax_table(ps_rel_pruned)[5,4] # let's say ASV5, so Micrococcales

# rename the NA's on family level as Order + Unknown
tax_table(ps_rel_pruned)[5,4] <- gsub("NA", "Micrococcales (Uknown)", tax_table(ps_rel_pruned)[5,4])

# You can do this on whatever taxonomic level you want, although if you try and do
# it on genus or species level it might be a lot depending on how many ASVs you
# choose to include

# Plot heat map using the Phyloseq package. Choose a transformation of the data,
# for example in log form. If you choose to order the taxa per family or whatever other label you choose 
# and samples per sample ID, it will order the heatmap the way you indicate. You can also let it set
# the order for you through some ordination form (CLR, bray-curtis etc.)

heatmap <- plot_heatmap(ps_rel_pruned, 
                        sample.label = "Sample_ID", taxa.label = "Family", trans = "log10",
                        low="#D7DBCF", high="#384E14", na.value="#D7DBCF", taxa.order = "Family", 
                        sample.order = "Sample_ID", stat = "unique", width = 0.2) +
  theme(axis.title.x = element_text(size=8, margin = margin(t=3)),
        axis.title.y = element_text(size=9, margin = margin(r=2)),
        axis.text.y  = element_text(size=6, color = "#3b3b3b", face = "italic"), 
        axis.text.x  = element_text(size=5, color = "#3b3b3b")) +
  facet_grid(~Sample_type, scales="free", space="free") + 
  theme(strip.text.x = element_text(size = 7, color = "white"), strip.background = element_rect(fill = "#42473a"))
heatmap$scales$scales[[2]]$name <- "Families" # change name of axis
heatmap$scales$scales[[1]]$name <- "Samples"
print(heatmap)

##### Core microbiome, using microbiome package ################################
# I vaguely follow this tutorial https://microbiome.github.io/tutorials/Core.html

# It's actually a good idea to make a new object that is normalized, this is the
# same as calculating the relative abundance
ps.rel <- microbiome::transform(ps_Lagilis, 'compositional')

# Re-name the taxonomic levels (you need to do this to edit the tax table)
# unfortunately the package uses different names so you need to change it
colnames(tax_table(ps.rel)) <-  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# This adds the best taxonomic assignment to the ASV; eg. ASV:Roseburia
ps_rel.f <- microbiome::add_besthit(ps.rel, sep = ":")

head(tax_table(ps_rel.f)) #check

#### Plot on genus level ####

# merge all ASVs belonging to the same genus
ps_rel_gen <- aggregate_taxa(ps_rel.f, "Genus")

# This will just set the x-axis, what increments you want the thresholds to be at
# for example, here my minimum abundance is 0.1%
detections <- round(10^seq(log10(0.001), log10(0.3), length = 15), 3)

taxa_names(ps_rel_gen) # get names

# Rename to Family (unknown), you will probably need to do this several times
taxa_names(ps_rel_gen) <- gsub("Bacteria_Actinobacteriota_Actinobacteria_Corynebacteriales_Corynebacteriaceae_NA",
                               "Corynebacteriaceae (Uknown)", taxa_names(ps_rel_gen))

# I set new color palette, chose your own per your tastes
colors2 = c("#262315", "#333200", "#565406", "#565406", "#8c8153", "#f2eee0")

# Now plot the core microbiota. In this case I have a minimum sample prevalence
# of 50% (0.5), i.e ASVs that are present in at least half of the samples. 
# You can change this per published data comparison for example.
# 0.7 is also a common threshold but quite strict

plot_core(ps_rel_gen, plot.type = "heatmap",
          detections = detections, min.prevalence = 0.5, colours = rev(colors2)) +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  theme_bw() + ylab("Genus") +
  theme(axis.text.x = element_text(angle=0),
        axis.text.y = element_text(face = "italic"))

############### Venn diagram of core microbiome ################################

# Chose a variable you want to compare, in my case sample type
table(meta(ps_rel.f)$Sample_type, useNA = "always")

# I just re-named a new object but you don't have to
ps_rel_venn <- ps_rel.f

# Make a vector which has all the sample types in it as a list
sample_types <- unique(as.character(meta(ps_rel_venn)$Sample_type))
print(sample_types)

list_core <- c() # an empty object to store information

for (n in sample_types){ # for each variable n in sample_types
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps_rel_venn, Sample_type == n) # Choose sample from sample type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 50% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(venn(list_core),
     fills = color_list[c(6,9,11)])
#main = list(label = "Shared and unique core ASVs (min. 50% prevalence)", cex = 1.4, vjust = 1, font = 2))

################################################################################
############################### End of pipeline ################################
################################################################################

# There are of course many other types of analysis that you can do, and this
# script is not comprehensive. Feel free to edit it and add in new things.