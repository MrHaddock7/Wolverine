# Libraries / packages

library(phyloseq)

# Loading data

source("config.R")

path_ps <- paste(home, "data/ps_1126.RData", sep="")
save_path <- file.path(home, "data")

load(path_ps)
ps <- ps_1126

######################## Make rarefaction curve ################################

# You can plot how many reads you need before most of the ASVs in a sample are
# present, this guides how many reads to rarefy with

# I chose to limit number of reads because the plot became quite illegible
# however adapt this to your data set
ps_rarecurve = prune_samples(sample_sums(ps)<=10000, ps)

tab <- otu_table(ps_rarecurve)
class(tab) <- "matrix" # as.matrix() will do nothing, you get a warning here, but this is what we need to have
tab2 <- t(tab) # transpose observations to rows

# Plot rarefaction curve, see where it levels out. It seems that my data set seems to
# level out at around approx. 1000-2000 reads, but of course this will differ
rare <- rarecurve(tab, step=1000, lwd=2, ylab="OTU",  label=F)

########################### Rarefy reads #######################################

# This takes x random reads from all samples (lowest number chosen) (which is
# why we set seed), to control for different number of reads for samples

# Alpha diversity calculations are sensitive to number of reads (makes sense),
# so then creating an even (random) read depth reduces biases

# Call it something that makes sense, in my case all my rarefied data is called
# ps_rare in some form, all my normalized data is called ps_rel (relative abundance)

set.seed(7444) # chose random number
# Set it to the lowest amount of reads you have, depends on your data set,
# also can be guided by rarefaction curve

# Checking lowest amount of reads = 7444, 12910 and 188325, we'll run it with 
# 188325
min(sample_sums(ps))

sort(sample_sums(ps))[1:10]   # lowest 10

# Removing the two low read samples ( P34104_331, P34104_333)

ps_filtered <- prune_samples(sample_sums(ps) > 12910, ps)
min(sample_sums(ps_filtered))

# Rarefying

ps_rare = rarefy_even_depth(ps, sample.size = 188325)

# Double checking
ntaxa(ps)
ntaxa(ps_rare)

# Now normalize by calculating the relative abundance of all the samples
ps_rel = transform_sample_counts(ps_filtered, function(x) x / sum(x) )

# Saving files, remember to change date
save(ps_rare, file = file.path(save_path, "ps_rare_1209.RData"))
save(ps_rel, file = file.path(save_path, "ps_rel_1209.RData"))
save(ps_filtered, file = file.path(save_path, "ps_filtered_1209.RData"))


# ps_filtered = bad samples removed
# ps_rare = for alpha diversity
# ps_rel = for barplots and composition

