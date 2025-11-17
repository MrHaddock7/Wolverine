library("dada2")
load("savepoint2.RData")

taxa <- assignTaxonomy(seqtab.nochim, "./silva/silva_trainset.fa.gz", multithreading=TRUE)
taxa <- addSpecies(taxa, "./silva/silva_v138.2_assignSpecies.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save(taxa, seqtab.nochim, file="taxa.RData")
