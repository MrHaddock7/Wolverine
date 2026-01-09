## This script is used to create a decontaminated phyloseq object using
## the DADA2 ASV output as an input. Contaminants, low-quality samples
## were filtered out for further downstream analysis.

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

install.packages("microbiomeutilities")
library(microbiomeutilities)


# Load packages, check which version you have 

library(dada2); packageVersion("dada2")
library(ShortRead)
library(ggplot2)
library(MASS)

library(phyloseq)
library(Biostrings)
library(decontam)
library(microbiome)
library(vegan)

# Some extra packages to use for plotting, should be easy to find and directly install
library(ggplot2)
library(grid)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(ggsignif)

### 1. Load metadata ----------------------------------------------------------

source("config.R")

path_metadata <- paste(home, "data/metadata.csv", sep="")

metadata <- read.csv(
  path_metadata,
  sep = ";",
  stringsAsFactors = FALSE
)

### 2. Extracting and cleaning sample names so they match the metadata identifiers

# Set rownames to NGI.ID
rownames(metadata) <- metadata$NGI.ID

# Extract sample names from seqtab.nochim 
samples.out <- rownames(seqtab.nochim)

# Remove trailing underscore
samples.clean <- sub("_$", "", samples.out)

# Remove Illumina sample index (_S###)
samples.clean <- sub("_S[0-9]+$", "", samples.clean)

# Now assign cleaned names to seqtab.nochim
rownames(seqtab.nochim) <- samples.clean

# Filter metadata to only samples that have ASV data 
metadata.filtered <- metadata[samples.clean, ]

# Ensure rownames match
rownames(metadata.filtered) <- samples.clean

### 3. Build phyloseq object --------------------------------------------------

ps_raw <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(metadata.filtered),
  tax_table(taxa)
)

### 4. Convert sequences into ASV names ---------------------------------------

dna <- Biostrings::DNAStringSet(taxa_names(ps_raw))
names(dna) <- taxa_names(ps_raw)
ps_raw <- merge_phyloseq(ps_raw, dna)

taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))

### 5. Clean up workspace (optional) ------------------------------------------

rm(
  dadaFs, dadaRs, errF, errR, mergers, out_trim, seqtab, seqtab2, track2,
  filtFs, filtRs, fnFs, fnRs, FWD, FWD.orients, path, REV, REV.orients
)

ps_raw

### 6. Identify and remove unwanted taxa ------------------------------------------

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

### 7. Prepare metadata for detecting contaminants ------------------------------------

# Check metadata
head(sample_data(ps_raw))

# Convert sample_data to data.frame for plotting
df <- as.data.frame(sample_data(ps_raw))
df$LibrarySize <- sample_sums(ps_raw)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Converting sample data to data frame and converting User_ID column to UTF-8
# Ensure User_ID is character
df$User_ID <- as.character(df$User_ID)

# Convert to UTF-8 to avoid translation warnings
df$User_ID <- iconv(df$User_ID, from = "", to = "UTF-8")

# Create Sample_or_Control column based on User_ID
df$Sample_or_Control <- ifelse(grepl("^Control_", df$User_ID),
                               "Control",
                               "Sample")

### 8. Identifying and visualizing contaminants -------------------------------------

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

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_raw, function(abund) 1*(abund > 0))

# Split into negative controls and true samples
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

### 9. Removing contaminants ------------------------------------------

# Prune contaminants
ps_noncontam <- prune_taxa(!contamdf.prev$contaminant, ps_raw)

# Remove the 'is.neg' column from sample_data
sample_data(ps_noncontam)$is.neg <- NULL

# Check final sample_data
sample_data(ps_noncontam)

# Clean up workspace
rm(df, df.pa, ps.pa, ps.pa.neg, ps.pa.pos, contamdf.prev)

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
ps_no_ctrl

# Check to see how much was removed
n_raw  <- ntaxa(ps_raw)
n_clean <- ntaxa(ps_no_ctrl)
n_raw
n_clean
removed <- n_raw - n_clean
removed
removed / n_raw * 100

# Choose a path and filename

path_save <- paste(home, "data/ps_no_ctrl_clean.RData", sep="")

save(ps_no_ctrl, file = path_save)
