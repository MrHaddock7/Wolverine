library(ggpicrust2)
library(phyloseq)
library(dplyr)




metadata <- as.data.frame(sample_data(ps_no_ctrl))
otu <- as.data.frame(otu_table(ps_no_ctrl))

pathway_abundance <- example$KO_pathway

# Perform DAA
daa_result <- perform_daa(
  input.data = pathway_abundance,
  sample.metadata = metadata,
  group = "Group"
)

# Visualize DAA results
visualize_daa(daa_result)

