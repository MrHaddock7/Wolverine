
# Script for providing a pie chart for relative abundance of specific sample
# the script does these things:
# 1. Normalises the phyloseq object (ps_rel) to relative abundance (compositional)
# 2. For each specific individual:
#     - Extracts individual's data
#     - Generates two pie charts (phylum & family)
#     - Collapses low abundance phyla into "Other"
#     - Optionally saves each pie as a .png
# 3. Returns all generated pie charts as a named list

# Usage example:
#   make_pie_plot(ps_rel, sample_ids = "P34104_331")
#   make_pie_plot(ps_rel, sample_ids = c("P34104_331", "P34104_332"), save_plot = TRUE)

# Testing and installing packages"
install.packages("patchwork")
install.packages("rlang")


make_pie_plot <- function(ps_rel,
                          sample_ids,
                          top_n = 5,
                          save_plot = FALSE,
                          out_dir = ".") {
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(microbiome)
  library(patchwork)
  
  # 1. Normalising
  ps_norm <- microbiome::transform(ps_rel, "compositional")
  df <- psmelt(ps_norm)
  
  # Check that the requested samples actually exist 
  missing <- setdiff(sample_ids, unique(df$NGI.ID))
  if (length(missing) > 0) {
    stop(paste("Sample IDs not found:", paste(missing, collapse = ",")))
  }
  
  # Storing all pie chart objects
  output_list <- list()
  
  # -- Function for one Pie -- #
  make_single_pie <- function(df_sub, tax_level, sample_id) {
    
    if (!tax_level %in% colnames(df_sub)) {
      stop(paste("Taxonomy level:", tax_level, "not found"))
    }
    
    # Identify top taxa
    taxa_summary <- aggregate(Abundance ~ ., 
                              data = df_sub[c(tax_level, "Abundance")], 
                              sum)
    taxa_summary <- taxa_summary[order(-taxa_summary$Abundance), ]
    top_taxa <- as.character(taxa_summary[[tax_level]][1:min(nrow(taxa_summary), top_n)])
    
    # Collapse rare taxa into "Other"
    df_sub[[tax_level]] <- as.character(df_sub[[tax_level]])
    df_sub[[tax_level]][!(df_sub[[tax_level]] %in% top_taxa)] <- "Other"
    df_sub[[tax_level]][is.na(df_sub[[tax_level]])] <- "Unassigned"
    
    # Summarise collapsed abundances
    df_clean <- aggregate(Abundance ~ ., 
                          data = df_sub[c(tax_level, "Abundance")], 
                          sum)
    df_clean$Percent <- df_clean$Abundance * 100
    
    # Generate colors for top taxa dynamically
    top_colors <- c("#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#E69F00")
    tax_levels <- unique(df_clean[[tax_level]])
    tax_colors <- setNames(rep("#999999", length(tax_levels)), tax_levels)  # default grey
    
    # Assign colors to top taxa
    for (i in seq_along(top_taxa)) {
      if (i <= length(top_colors)) {
        tax_colors[top_taxa[i]] <- top_colors[i]
      }
    }
    
    # Assign specific colors for Other and Unassigned
    if ("Other" %in% tax_levels) tax_colors["Other"] <- "#999999"
    if ("Unassigned" %in% tax_levels) tax_colors["Unassigned"] <- "#CC79A7"
    
    # Plot
    pie <- ggplot(df_clean, aes(x = "", y = Abundance, fill = !!sym(tax_level))) +
      geom_col(width = 1, color = "white") +
      coord_polar(theta = "y") +
      theme_void() +
      ggtitle(paste("Relative Abundance (", tax_level, "): ", sample_id, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            legend.title = element_blank()) +
      scale_fill_manual(values = tax_colors)
    
    
    return(pie)
  }
  
  # -- Loop function for samples -- #
  for (sample_id in sample_ids) {
    
    df_sub <- df %>% filter(NGI.ID == sample_id)
    
    # Creating pies
    pie_phylum <- make_single_pie(df_sub, "Phylum", sample_id)
    pie_family <- make_single_pie(df_sub, "Family", sample_id)
    
    # Combining side-by-side
    combined <- pie_phylum + pie_family + plot_layout(ncol = 2)
    
    # Saving combined plot
    if (save_plot) {
      out_combined <- file.path(out_dir, paste0("pie_", sample_id, "_combined.png"))
      ggsave(out_combined, combined, width = 10, height = 6)
    }
    
    # Saving plots if requested
    if (save_plot) {
      out1 <- file.path(out_dir, paste0("pie_", sample_id, "_phylum.png"))
      out2 <- file.path(out_dir, paste0("pie_", sample_id, "_family.png"))
      
      ggsave(out1, pie_phylum, width = 6, height = 6)
      ggsave(out2, pie_family, width = 6, height = 6)
      
      message("Saved: ", out1)
      message("Saved: ", out2)
    }
    
    # Store in the output list 
    output_list[[paste0(sample_id, "_Phylum")]] <- pie_phylum
    output_list[[paste0(sample_id, "_Family")]] <- pie_family
    output_list[[paste0(sample_id, "_Combined")]] <- combined
  }
  
  return(output_list)
}


pies <- make_pie_plot(ps_rel, sample_ids = "P34104_353")
pies$P34104_353_Phylum
pies$P34104_353_Family
pies$P34104_353_Combined

pies <- make_pie_plot(ps_rel, 
                      sample_ids = c("P34104_353", "P34104_354"))
pies[["P34104_353_Phylum"]]     # Phylum pie for sample P34104_353
pies[["P34104_353_Family"]]     # Family pie for sample P34104_353
pies[["P34104_353_Combined"]]   # Combined plot for sample P34104_353

pies[["P34104_354_Phylum"]]     # Phylum pie for sample P34104_354
pies[["P34104_354_Family"]]     # Family pie for sample P34104_354
pies[["P34104_354_Combined"]]   # Combined plot for sample P34104_354


