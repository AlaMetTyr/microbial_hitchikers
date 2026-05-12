library(phyloseq)

# Subset taxa where species is muscaria
ps_amanita <- subset_taxa(ps_final, Genus == "g__Amanita")

# Extract OTU/ASV table
otu_ama <- otu_table(ps_amanita)

# Convert to presence/absence
presence <- otu_ama > 0

# Get sample names where present
samples_positive_ama <- sample_names(ps_amanita)[colSums(presence) > 0]

samples_positive_ama


###try with asv counts
target_asvs <- c("ASV_272", "ASV_564", "ASV_2241", "ASV_2352")

ps_subset <- prune_taxa(target_asvs, ps)

otu <- otu_table(ps_subset)

samples_positive <- sample_names(ps_subset)[colSums(otu > 0) > 0]

samples_positive
