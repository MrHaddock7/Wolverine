library("dada2")
load("read_savepoint2.RData")

taxa <- assignTaxonomy(seqtab.nochim, "~/silva_trainset.fa.gz", multithreading=FALSE)
taxa <- addSpecies(taxa, "/mnt/c/Users/Lovisa/Documents/Skola/Tillampad_Bioinformatik/silva_v138.2_assignSpecies.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save(taxa, file="taxa.RData")
