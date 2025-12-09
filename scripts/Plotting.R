# Libraries / packages

library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)  # for stat_compare_means / stat_kruskal_test


################################################################################
######################### Plotting and data-exploration ########################
################################################################################

# Loading data

source("config.R")

path_ps <- paste(home, "data/ps_1126.RData", sep="")

load(path_ps)
ps <- ps_1126

# Set color-palette for plotting

color_list <- c(
  "Bacillota"      = "#FF9999",  # light red
  "Pseudomonadota" = "#99CCFF",  # light blue
  "Actinomycetota" = "#66CC99",  # teal
  "Bacteroidota"   = "#FFCC66",  # light orange
  "Chloroflexota"  = "#CC99CC",  # pastel purple
  "Other"          = "#CCCCCC"   # grey
)

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

set.seed(300) # chose random number
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

# Saving ps_rare as file
save(ps_rare, file = "ps_rare.RData")

ntaxa(ps)
ntaxa(ps_rare)

# Now normalize by calculating the relative abundance of all the samples
ps_rel = transform_sample_counts(ps_filtered, function(x) x / sum(x) )

# ps_filtered = bad samples removed
# ps_rare = for alpha diversity
# ps_rel = for barplots and composition

# I would save the work space at this stage and call it "Phyloseq_decontam" or
# similar, so you can always re-load from this stage

################################################################################
######## Relative abundance box-plots ###########################################
################################################################################

# This is for the bar-plots showing relative abundance of all samples #
# I'll give example on phylum level, but it's then easy to change the taxonomic level 
# and re-do it on family or genus level for example

# Normalize
ps_norm <- microbiome::transform(ps_rel, "compositional")
barplot_phylum <- psmelt(ps_norm)

# Replacing Sample_Type with NGI.ID
# Turn sample_ID into factor, since they are numbers but don't want them to behave as numbers
# or whatever you called that column
# also turn whatever variable you're interested in comparing (for example location)
# into a factor! For me it was sample type
barplot_phylum$NGI.ID <- as.factor(barplot_phylum$NGI.ID)  
barplot_phylum$Zoo <- as.factor(barplot_phylum$Zoo)

# Top phyla
taxa_summary <- aggregate(Abundance ~ Phylum, data = barplot_phylum, sum)
taxa_summary <- taxa_summary[order(-taxa_summary$Abundance), ]
top_phyla <- as.character(taxa_summary$Phylum[1:5])
top_phyla

# For-loop and if then function that checks to see if any of the names in the list "top" is in the 
# barplot_phylum$Phylum column, leaves unchanged if it is, changes it to "Other" if it isn't
# Re-name all uncommon phyla as "Other"
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in seq_along(barplot_phylum$Phylum)) {
  if (barplot_phylum$Phylum[i] %!in% top_phyla) {
    barplot_phylum$Phylum[i] <- "Other"
  }
}

# Re-name NA:s as "Other"
barplot_phylum$Phylum[is.na(barplot_phylum$Phylum)] <- "Unassigned"
barplot_phylum$Phylum <- factor(barplot_phylum$Phylum, levels= c(top_phyla, "Other"))

# Change order so "Other" is last, change the rest to whichever order you prefer
barplot_phylum$Phylum <- factor(barplot_phylum$Phylum, levels= c("Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota","Chloroflexota", "Other"))
levels(barplot_phylum$Phylum) # check

# Plot per individual sample sorted by sample type (or other variable)
# In my case I have "Sample_ID", which is the number for all my samples, 
# but your data set might have this organized differently
relab_plot <- ggplot(barplot_phylum, aes(x = NGI.ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") +
  facet_grid(~Zoo, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
relab_plot
# I am also using something called "facet_grid", that separates out my variable of
# interest (sample type) into different grids.

# Saving plot
ggsave("relative_abundance_plot.pdf", plot = relab_plot, width = 12, height = 6, units = "in")

########################## Merged barplot #######################################

# This might be of interest for you, if you want to merge the samples by 
# one variable, for example sex or location or in my case sample type

colnames(sample_data(ps_norm))

# Merging all samples with the same Zoo value 

sample_data(ps_norm)$Zoo <- as.factor(sample_data(ps_norm)$Zoo)
table(sample_data(ps_norm)$Zoo, useNA="ifany")

# Merge raw counts by Zoo
ps_merged_phylum <- merge_samples(ps, "Zoo")  

# Normalize after merging
merged_phylum_rel <- transform_sample_counts(ps_merged_phylum, function(x) x / sum(x))

# Melt for plotting
barplot_phylum_merged <- psmelt(merged_phylum_rel)

# Merge rare phyla into "Other"
barplot_phylum_merged$Phylum[barplot_phylum_merged$Phylum==""] <- NA
'%!in%' <- function(x,y)!('%in%'(x,y))
for(i in seq_along(barplot_phylum_merged$Phylum)) {
  if(barplot_phylum_merged$Phylum[i] %!in% top_phyla) {
    barplot_phylum_merged$Phylum[i] <- "Other"
  }
}
barplot_phylum_merged$Phylum[is.na(barplot_phylum_merged$Phylum)] <- "Other"
barplot_phylum_merged$Phylum <- factor(barplot_phylum_merged$Phylum, 
                                       levels=c("Bacillota","Pseudomonadota","Actinomycetota",
                                                "Bacteroidota","Chloroflexota","Other"))

# Plot using the Sample column already in barplot_phylum_merged
merged_relab_plot_zoo <- ggplot(barplot_phylum_merged, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  labs(x = "Zoo", y = "Relative Abundance") +
  scale_y_continuous(limits = c(0, 1))

merged_relab_plot_zoo

# Saving plot
ggsave("relative_abundance_merged_by_zoo.pdf", plot = merged_relab_plot_zoo, width = 12, height = 6, units = "in")

################################################################################
###### Plot different diversity indices ########################################
################################################################################

# Make box-plots of various diversity indices, here I compare sample type but 
# change it to whatever else you want. The plots are somewhat ugly so if you
# want to use it you may need to edit it more.

Zoo_div_faceted <- plot_richness(ps_rare, x = "Zoo", measures = c("Shannon", "Simpson"), color = "Zoo") +
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
ggsave("shannon_simpson_zoo.pdf", plot = Zoo_div_faceted, width = 12, height = 6, units = "in")

################################################################################
########## Box-plots for ASV richness ##########################################
################################################################################

# Since we're using alpha diversity it's good to use the rarefied object

# Zoos
Zoo_boxplot <- plot_richness(ps_rare, x = "Zoo", measures = "Observed")

ggplot(Zoo_boxplot$data, aes(x = Zoo, y = value, fill = Zoo)) +
  geom_boxplot(alpha = 1, color = "black", outlier.color = "white", outlier.size = 0, 
               lwd=0.36, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.25, color = "black", pch = 21) +
  ylab("Alpha Diversity") + xlab("Zoo") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "grey29"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.key.size = unit(1,"cm")) +
  theme(legend.position = "none") +
  stat_kruskal_test(group.by = "x.var", label = "p.signif", label.y = 220) +
  # Note here that I add in a statistical test to indicate significance
  # It depends on your data which test you should use
  ggtitle("Difference in ASV richness between Zoos")

#wilcox pairwise

# Sex
Sex_boxplot <- plot_richness(ps_rare, x = "Sex", measures = "Observed")

sample_data(ps_rare)$Sex <- factor(sample_data(ps_rare)$Sex,
                                   levels = c(0,1),
                                   labels = c("M", "F"))


ggplot(Sex_boxplot$data, aes(x = Sex, y = value, fill = Sex)) +
  geom_boxplot(alpha = 1, color = "black", outlier.color = "white", outlier.size = 0, 
               lwd=0.36, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.25, color = "black", pch = 21) +  
  ylab("Alpha Diversity (Observed ASVs)") + xlab("Sex") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "grey29"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_fill_manual(values = color_list[c(6,4,8)]) +
  scale_color_manual(values = color_list) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.key.size = unit(1,"cm")) +
  theme(legend.position = "none") +
  stat_kruskal_test(group.by = "x.var", label = "p.signif", label.y = 220) +
  # Note here that I add in a statistical test to indicate significance
  # It depends on your data which test you should use
  ggtitle("Difference in ASV richness between Sex")

# Diet
Diet_boxplot <- plot_richness(ps_rare, x = "Diet", measures = "Observed")

ggplot(Diet_boxplot$data, aes(x = Diet, y = value, fill = Diet)) +
  geom_boxplot(alpha = 1, color = "black", outlier.color = "white", outlier.size = 0, 
               lwd=0.36, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.25, color = "black", pch = 21) +  
  ylab("Alpha Diversity (Observed ASVs)") + xlab("Diet") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "grey29"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_fill_manual(values = color_list[c(6,4,8)]) +
  scale_color_manual(values = color_list) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.key.size = unit(1,"cm")) +
  theme(legend.position = "none") +
  stat_kruskal_test(group.by = "x.var", label = "p.signif", label.y = 220) +
  # Note here that I add in a statistical test to indicate significance
  # It depends on your data which test you should use
  ggtitle("Difference in ASV richness between Diet")

# Example on how you might plot two different box-plots in the same plot
# using ggarrange. You can't use par(mfrow = c(2,2)) with ggplot 
ggarrange(plot1_name, plot2_name, nrow = 1, labels = c("A", "B"))

################################################################################
################################# PCoA/PCA plot ################################
################################################################################

# Using a PCoA (principal coordinate analysis) or PCA (principal component analysis)
# is a nice way to visualize the ASV composition of your samples, i.e the beta diversity

# I ended up using Aitchinson's distance but I added examples for both
# Bray-Curtis ordination and Aitchinson's distance.
# It might be good to try different methods.

# Turn your sample variables into factors
# Important that you use the normalized data set:
sample_data(ps_norm)$Zoo <- as.factor(sample_data(ps_norm)$Zoo)
sample_data(ps_norm)$Sex <- as.factor(sample_data(ps_norm)$Sex)


################################################################################
# Elin's PcOA plot example

# change here the colors to what you like or want to use. Use as many colors as you have categories. 
col_gir[2] <- "#c45c1a"
col_gir[1] <- "#F2C188"
names(col_gir) <- c("Northern", "Reticulated", "Masai")

ord_ps = ordinate(psF.rel, "PCoA", "bray")
g <- plot_ordination(psF.rel, ord_ps, justDF=T)
head(g)

# extract % variation
tempx <- as.character(round(ord_ps$values$Relative_eig[1],3) *100)
tempy <- as.character(round(ord_ps$values$Relative_eig[2],3) *100)
set.seed(386)

# pcoa plot Giraffe Species
ggplot(g, aes(x=Axis.1, y=Axis.2, color=Species, fill=Species)) + theme_bw() + 
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, level=0.9)  +  # confidence interval
  geom_point(color="black", shape=21, size = 2.5, stroke=0.3) +
  xlab(paste0("PC1 (",temps,"%)")) + ylab(paste0("PC2 (",tempy,"%)")) +
  theme(axis.title.x = element_text(size=13, margin = margin(t=3)),
        axis.title.y = element_text(size=13, margin = margin(r=2)),
        axis.text.y  = element_text(size=11, color = "black"), 
        axis.text.x  = element_text(size=11, color = "black"),
        strip.text   = element_text(size=13)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 13), 
        legend.text  = element_text(size = 13),
        legend.key.size = unit(0.6, "cm"),
        legend.position = "right") +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values  = color_list)

################################################################################
# PCA plot example

# This transforms the relative abundances using "clr" (centered log ratio) transformation
ps_clr <- microbiome::transform(ps_rel, "clr")

# Ordinate RDA (redundancy analysis), this is for the distances for the PCA
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

# calculate the relative eigenvalues for pc1 and pc2
# $CA is the object in ord_clr that includes the eigenvalues for each PC
sum1 <- sum(ord_clr$CA$eig)

tempx <- as.character(round((ord_clr$CA$eig[1]/sum1)*100))
tempy <- as.character(round((ord_clr$CA$eig[2]/sum1)*100))

PCA_aitch <- plot_ordination(ps_rel, ord_clr, justDF=T)

ggplot(PCA_aitch, aes(x=PC1, y=PC2, fill = Sample_type, color = Sample_type, group = Sample_type)) + 
  theme_bw() + 
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, level=0.95) +  # confidence interval
  geom_point(aes(shape = Sex), size = 2.5, stroke=0.3) +  # move shape aesthetic here
  labs(x=(paste0("PC1 (",tempx,"%)")), y=(paste0("PC2 (",tempy,"%)"))) +
  theme(axis.title.x = element_text(size=13, margin = margin(t=3)),
        axis.title.y = element_text(size=13, margin = margin(r=2)),
        axis.text.y  = element_text(size=11, color = "black"), 
        axis.text.x  = element_text(size=11, color = "black"),
        strip.text   = element_text(size=13)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 13), 
        legend.text  = element_text(size = 13),
        legend.key.size = unit(0.6, "cm"),
        legend.position = "right") +
  scale_color_manual(values = color_list[c(9,5,1)]) +
  scale_fill_manual(values = color_list[c(9,5,2)], guide = "none") +
  scale_shape_manual(values = c(21,22)) +  # specify shape values
  guides(color = guide_legend("Sample type")) # this is to make the legend colors work

################################################################################
############## Statistical tests on alpha diversity ############################
################################################################################

# These are tests to statistically compare the alpha diversity of different variables in the
# data. You might need to do different tests depending on your data. 

# It is probably best to use non-parametric tests with microbiome data, although you 
# can test for normality by plotting Q-Q plots, doing a quick histogram or a shapiro test. 
# A p-value > 0.05 means it's normally distributed

# hist(variable)
# qqplot(variable)
# qqline(variable)
# shapiro.test(variable)

# Get ASV richness data, you can use the "estimate_richness" function
# you could also use a different alpha diversity measure than observed richness
rich = estimate_richness(ps_rare)

# Here you can go back and use your metadata, as it will be nicely organized
# and contain all the information you'll need
# you will need to remove all the actual samples you removed before, e.g.
# your control samples and anything else you removed. 

alpha_div = subset(metadata, Sample_or_control == "Sample")

# check to see that row numbers are the same
nrow(alpha_div)
nrow(rich)

# I added an extra column with the observed richness values, you will
# need to check how many columns you have 
ncol(alpha_div)
alpha_div[,15] <- rich$Observed
names(alpha_div)[15] <- paste("Richness")

# Here I have an example of a Kruskal-Wallice test, which is the non-parametric
# one-way ANOVA, testing for two or more groups on a continuous variable

# Here I test sample types against each other, but of course change the variable
# per your data set
kw1 <- kruskal.test(Richness ~ Sample_type, data = alpha_div)

# Plot distribution, you can transform data if it's not normal, or change
# the distribution in the glm later on.
hist(alpha_div$Richness)

# And here I test a general linear model. Change the family to whatever distribution
# you have. It could for example be gaussian, poisson etc.
glm_1 <- glm(Richness ~ Sex + Temp_group, data = alpha_div, family = gaussian)
summary.glm(glm_1)

################################################################################
########## Adonis tests on ASV composition #####################################
################################################################################

# The Adonis test is included in the microbiome package, and is a version of a 
# PERMANOVA (Permutational analysis of variance). It is essentially a test for
# difference in beta diversity between sample variables. 

# Begin by using your normalized data (since we're testing beta diversity),
# and all your variables need to be as factors, so transform them as so
sample_data(ps_rel)$Sex <- as.factor(sample_data(ps_rel)$Sex)

############################## Using CLR #######################################

# centered log-ratio (CLR) transformation of counts
ps_clr <- microbiome::transform(ps_rel, "clr")

#Generate distance matrix, you can do this in various ways. 
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean")

#adonis test, on sample type, sex, location etc, whatever variable you want
# you will get a P-value, and an R2 value (how much variation is explained)
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Sex, permutations = 999)

######################## Using bray-curtis #####################################

#Generate distance matrix, you can do this in various ways. 
bray_dist_matrix <- phyloseq::distance(ps_rel, method = "bray")

#adonis test, on sample type, sex, location etc, whatever variable you want
# you will get a P-value, and an R2 value (how much variation is explained)
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps_rel)$Sex, permutations = 999)

# You can also run it on two variables, if so add a "+" between them
# if you use * you also test the interaction, which is more complicated to interpret

################################################################################
############################ Core microbiome ################################### 
################################################################################

# https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_heatmap
# https://joey711.github.io/phyloseq/plot_heatmap-examples.html

# The core microbiome can be defined in different ways. It is also somewhat
# ambiguous what the biological meaning actually is. It can however be useful
# if you want to compare your data directly to someone else's study, so I have
# included it here.

# First let us just create a regular heatmap plot of the most abundant taxa. You
# may want this to show the composition, so it's not actually the core microbiota
# yet. 

# Filter out 50 most prevalent ASVs (or however many you want)
ps_rel_pruned <- prune_taxa(names(sort(taxa_sums(ps_rel),TRUE)[1:50]), ps_rel)

tax_table(ps_rel_pruned) # get tax table for the data set

# Depending on if you have a lot of NA:s in your plot, you may want to re-name
# them to something nicer, here I have an example:

print(tax_table(ps_rel_pruned)[,5]) # this shows all top ASVs at the genus level 

# if any of the ASVs have NA:s note them down

# find the name of the order the unknown family/genus belongs to 
tax_table(ps_rel_pruned)[5,4] # let's say ASV5, so Micrococcales

# rename the NA's on family level as Order + Unknown
tax_table(ps_rel_pruned)[5,4] <- gsub("NA", "Micrococcales (Uknown)", tax_table(ps_rel_pruned)[5,4])

# You can do this on whatever taxonomic level you want, although if you try and do
# it on genus or species level it might be a lot depending on how many ASVs you
# choose to include

# Plot heat map using the Phyloseq package. Choose a transformation of the data,
# for example in log form. If you choose to order the taxa per family or whatever other label you choose 
# and samples per sample ID, it will order the heatmap the way you indicate. You can also let it set
# the order for you through some ordination form (CLR, bray-curtis etc.)

heatmap <- plot_heatmap(ps_rel_pruned, 
                        sample.label = "Sample_ID", taxa.label = "Family", trans = "log10",
                        low="#D7DBCF", high="#384E14", na.value="#D7DBCF", taxa.order = "Family", 
                        sample.order = "Sample_ID", stat = "unique", width = 0.2) +
  theme(axis.title.x = element_text(size=8, margin = margin(t=3)),
        axis.title.y = element_text(size=9, margin = margin(r=2)),
        axis.text.y  = element_text(size=6, color = "#3b3b3b", face = "italic"), 
        axis.text.x  = element_text(size=5, color = "#3b3b3b")) +
  facet_grid(~Sample_type, scales="free", space="free") + 
  theme(strip.text.x = element_text(size = 7, color = "white"), strip.background = element_rect(fill = "#42473a"))
heatmap$scales$scales[[2]]$name <- "Families" # change name of axis
heatmap$scales$scales[[1]]$name <- "Samples"
print(heatmap)

##### Core microbiome, using microbiome package ################################
# I vaguely follow this tutorial https://microbiome.github.io/tutorials/Core.html

# It's actually a good idea to make a new object that is normalized, this is the
# same as calculating the relative abundance
ps.rel <- microbiome::transform(ps_Lagilis, 'compositional')

# Re-name the taxonomic levels (you need to do this to edit the tax table)
# unfortunately the package uses different names so you need to change it
colnames(tax_table(ps.rel)) <-  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# This adds the best taxonomic assignment to the ASV; eg. ASV:Roseburia
ps_rel.f <- microbiome::add_besthit(ps.rel, sep = ":")

head(tax_table(ps_rel.f)) #check

#### Plot on genus level ####

# merge all ASVs belonging to the same genus
ps_rel_gen <- aggregate_taxa(ps_rel.f, "Genus")

# This will just set the x-axis, what increments you want the thresholds to be at
# for example, here my minimum abundance is 0.1%
detections <- round(10^seq(log10(0.001), log10(0.3), length = 15), 3)

taxa_names(ps_rel_gen) # get names

# Rename to Family (unknown), you will probably need to do this several times
taxa_names(ps_rel_gen) <- gsub("Bacteria_Actinobacteriota_Actinobacteria_Corynebacteriales_Corynebacteriaceae_NA",
                               "Corynebacteriaceae (Uknown)", taxa_names(ps_rel_gen))

# I set new color palette, chose your own per your tastes
colors2 = c("#262315", "#333200", "#565406", "#565406", "#8c8153", "#f2eee0")

# Now plot the core microbiota. In this case I have a minimum sample prevalence
# of 50% (0.5), i.e ASVs that are present in at least half of the samples. 
# You can change this per published data comparison for example.
# 0.7 is also a common threshold but quite strict

plot_core(ps_rel_gen, plot.type = "heatmap",
          detections = detections, min.prevalence = 0.5, colours = rev(colors2)) +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  theme_bw() + ylab("Genus") +
  theme(axis.text.x = element_text(angle=0),
        axis.text.y = element_text(face = "italic"))

############### Venn diagram of core microbiome ################################

# Chose a variable you want to compare, in my case sample type
table(meta(ps_rel.f)$Sample_type, useNA = "always")

# I just re-named a new object but you don't have to
ps_rel_venn <- ps_rel.f

# Make a vector which has all the sample types in it as a list
sample_types <- unique(as.character(meta(ps_rel_venn)$Sample_type))
print(sample_types)

list_core <- c() # an empty object to store information

for (n in sample_types){ # for each variable n in sample_types
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps_rel_venn, Sample_type == n) # Choose sample from sample type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 50% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(venn(list_core),
     fills = color_list[c(6,9,11)])
#main = list(label = "Shared and unique core ASVs (min. 50% prevalence)", cex = 1.4, vjust = 1, font = 2))

################################################################################
############################### End of pipeline ################################
################################################################################

# There are of course many other types of analysis that you can do, and this
# script is not comprehensive. Feel free to edit it and add in new things.