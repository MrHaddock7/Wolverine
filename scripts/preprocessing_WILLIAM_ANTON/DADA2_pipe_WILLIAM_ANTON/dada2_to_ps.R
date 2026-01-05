if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")

BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
BiocManager::install("Biostrings")
library(Biostrings); packageVersion("Biostrings")
BiocManager::install("ggplot2")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())


# 1. Read your metadata file
metadata <- read.csv("C:/Users/Lovisa/Documents/Skola/Tillampad_Bioinformatik/metadata.csv", stringsAsFactors = FALSE)

# 2. Ensure column name matches the rownames of seqtab.nochim
#    If it's called "NGI.ID" or similar, adjust here:
rownames(seqtab.nochim) <- sub("^([^_]+_[^_]+).*", "\\1", rownames(seqtab.nochim))

colnames(metadata)[colnames(metadata) == "NGI.ID"] <- "NGI_ID"


# 3. Set NGI_ID as rownames
rownames(metadata) <- metadata$NGI_ID

# 4. Subset metadata to only samples present in seqtab.nochim
samples.out <- rownames(seqtab.nochim)
samples.out
samdf <- metadata[samples.out, ]

samdf



ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps), ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
save(ps, file="ps.RData")
