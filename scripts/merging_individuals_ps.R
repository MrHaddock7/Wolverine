###Remove outliers
sample_names(ps_1126)
ps_filt <- prune_samples(!(sample_names(ps_1126) %in% c("P34104_331", "P34104_333")), ps_1126)
sample_names(ps_filt)


###Merge the samples
ps_merged <- merge_samples(ps_filt, "Individual")

meta <- data.frame(sample_data(ps_filt))

# Create one row per Individual
meta_merged <- meta[!duplicated(meta$Individual), ]

# Make Individual the rownames to match merged sample names
rownames(meta_merged) <- meta_merged$Individual

# Apply metadata to phyloseq object
sample_data(ps_merged) <- sample_data(meta_merged)


sample_names(ps_merged)
ps_merged_1201 <- ps_merged
save(ps_merged_1201, file = "C:/Users/Lovisa/Documents/Wolverine/data/ps_merged_1201.RData")
