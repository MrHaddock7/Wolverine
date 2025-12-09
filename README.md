# Wolverine

Data overview:

ps.RData = original phyloseq object
ps_clean_no_ctrl.RData = uses ps.RData to remove controls
ps_1126.RData = 
ps_rare.RData = rarefied ps_1126 for alpha diversity
ps_filtered.RData = removes low read samples _331 and _333
ps_rel.RData = normalised ps_filtered file by calculating the relative abundance 
              of all the samples


Scripts overview:

Piechart_script.R = to be used in automatic report generation, creates pie charts for showing
                    composition (Phylum/Family) in samples per zoo. 
                    
Rarefaction.R = rarefies data, results in ps_filtered, ps_rel and ps_rare files.

Relative_Abundance.R = Uses ps_rel to normalise data and generate relative abundance in sample
                        per zoo, also creates relative abundance plot for merged samples by zoo