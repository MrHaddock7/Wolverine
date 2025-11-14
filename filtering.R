library(dada2)
library(ShortRead)
library(Biostrings)

path <- "/home/haddock/private/järv/data/trimmed_primers"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "L001"), `[`, 1)

sample.names

# Ensure filtered folder exists
if(!dir.exists(file.path(path, "filtered"))) {
  dir.create(file.path(path, "filtered"))
}

#Inspect read quality profiles
#plotQualityProfile(fnFs)
#plotQualityProfile(fnRs)

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


length(filtFs) == length(unique(filtFs))

length(filtRs) == length(unique(filtRs))


table(filtFs)

table(filtRs)

length(unique(filtRs))
length(unique(filtFs))
unique(filtFs)
unique(filtRs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
    minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

save.image(file = "savepoint1.RData")
