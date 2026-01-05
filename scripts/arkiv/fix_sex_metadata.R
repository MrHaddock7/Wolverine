metadata <- read.csv("C:/Users/Lovisa/Documents/Wolverine/data/metadata.csv", sep = ";")
metadata$NGI.ID <- as.character(metadata$NGI.ID)
ps_sd <- data.frame(as(sample_data(ps_filtered), "data.frame"))
ps_sd$SampleID <- rownames(ps_sd) 
library(dplyr)


metadata
ps_sd2 <- ps_sd %>%
  left_join(metadata %>% select(NGI.ID, Sex), 
            by = c("SampleID" = "NGI.ID"))


rownames(ps_sd2) <- ps_sd2$SampleID

# 2. Remove the helper column
ps_sd2$SampleID <- NULL

# 3. Reorder rows to match phyloseq sample order
ps_sd2 <- ps_sd2[sample_names(ps_filtered), ]

# 4. Verify
identical(rownames(ps_sd2), sample_names(ps_filtered))

sample_data(ps_filtered) <- sample_data(ps_sd2)

table(sample_data(ps_filtered)$Sex)
identical(sample_names(ps_filtered), rownames(sample_data(ps_filtered)))
sample_data(ps_filtered)

save(ps_filtered, file = "C:/Users/Lovisa/Documents/Wolverine/data/ps_filtered_1209.RData")
