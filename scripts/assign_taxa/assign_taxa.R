library("dada2")
load("savepoint2.RData")

silva_train <- "/proj/rnaseq01/viroids/Viroid-transcriptome-mining/slutkurs/silva/silva_trainset.fa.gz"
silva_assign <- "/proj/rnaseq01/viroids/Viroid-transcriptome-mining/slutkurs/silva/silva_trainset.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithreading=TRUE)
taxa <- addSpecies(taxa, silva_assign)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save(taxa, file="taxa.RData")
