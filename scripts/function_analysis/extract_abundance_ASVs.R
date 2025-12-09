library(phyloseq)
library(Biostrings)

source("config.R")


path <- paste(home, "data/ps_1126.RData", sep="")
path

load(path)

seqs <- refseq(ps_1126)   # DNAStringSet med namn ASV1, ASV2, ...

# säkerställ att names(seqs) matchar taxa_names (ska de redan göra)
names(seqs)[1:5]
taxa_names(ps_1126)[1:5]

save_path <- paste(home, "data/ASV_abundance_sequences/rep_seqs_1126.fna", sep="")

# skriv fasta-fil till disk
writeXStringSet(seqs, filepath = save_path)



library(biomformat)

# Ta ut otutabellen som matrix
otu_mat <- as(otu_table(ps_1126), "matrix")

# PICRUSt2 vill ha taxa som rader: kontrollera
taxa_are_rows(ps_1126)
# Om FALSE, transponera:
if (!taxa_are_rows(ps_1126)) {
  otu_mat <- t(otu_mat)
}

# Skapa BIOM-objekt
otu_biom <- make_biom(data = otu_mat)

save_path_abundance <- paste(home, "data/ASV_abundance_sequences/table_1126.biom", sep="")
# Skriv till fil
write_biom(otu_biom, save_path_abundance)