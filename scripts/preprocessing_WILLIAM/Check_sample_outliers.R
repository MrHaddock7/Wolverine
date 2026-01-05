###This code was used to check for outliers among the samples, 
###We identified sample 331 and 333 as having low depth, ASV count 

library("phyloseq")
library("ggplot2")
sample_depths <- sample_sums(ps_no_ctrl)
sample_depths
bad_depth <- names(sample_depths[sample_depths < 15000])
bad_depth


###Check alpha diversity and observed ASV count per sample
library(phyloseq)
alpha <- estimate_richness(ps_1126, measures=c("Shannon", "Observed"))
alpha


ps.prop <- transform_sample_counts(ps_no_ctrl, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

p <- plot_ordination(ps.prop, ord.nmds.bray, color="Zoo", title="Bray NMDS")

p + geom_label(aes(label = NGI.ID), 
               size = 3, 
               alpha = 0.7)  + ggtitle(
                 paste0(
                   "Bray–Curtis NMDS (stress = ",
                   round(ord.nmds.bray$stress, 3),
                   ")"
                 )
               )
