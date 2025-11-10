library(dada2); packageVersion("dada2")

path <- '/Users/william/Library/CloudStorage/OneDrive-Uppsalauniversitet/Slutkurs/DADA2 tutorial/MiSeq_SOP'

list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)

# Forward. quality
plotQualityProfile(fnFs[1:2])

# Reverse quality
plotQualityProfile(fnRs[1:2])


# Filtering

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(240, 160), 
                     maxN = 0, 
                     maxEE = c(2, 2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

head(out)






