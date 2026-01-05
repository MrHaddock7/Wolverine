savepoint <- "/home/haddock/private/Wolverine/runs/run_trimming_taxa_EE2_4_17112025_57525403/savepoint2.RData"

library("dada2")
load(savepoint)

silva_train <- "/proj/rnaseq01/viroids/Viroid-transcriptome-mining/slutkurs/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_assign <- "/proj/rnaseq01/viroids/Viroid-transcriptome-mining/slutkurs/silva/silva_v138.2_assignSpecies.fa.gz"

x <- "training"
x
taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithread=TRUE)

x <- "assigning"
x
taxa <- addSpecies(taxa, silva_assign)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save(taxa, file="/home/haddock/private/Wolverine/runs/run_trimming_taxa_EE2_4_17112025_57525403/taxa.RData")

