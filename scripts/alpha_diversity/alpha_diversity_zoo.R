# Libraries

library(phyloseq)
library(ggplot2)
library(dplyr)

# Loading data

source("config.R")

path_ps  <- file.path(home, "data", "ps_rare_1209.RData")
save_path <- file.path(home, "plots_and_results")

load(path_ps)
ps <- ps_rel

################################################################################
###### Plot different diversity indices ########################################
################################################################################

# Make box-plots of various diversity indices, here I compare sample type but 
# change it to whatever else you want. The plots are somewhat ugly so if you
# want to use it you may need to edit it more.

Zoo_div_faceted <- plot_richness(ps, x = "Zoo", measures = c("Shannon", "Simpson"), color = "Zoo") +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_wrap(~variable, scales = "free_y")  # facet by diversity index

Zoo_div_faceted

# Saving plot
ggsave(
  filename = file.path(save_path, "alpha_diversity_shannon_simpson_zoo_1209.png"),
  plot = Zoo_div_faceted,
  width = 12,
  height = 6,
  units = "in"
)
