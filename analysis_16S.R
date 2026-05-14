library(tidyverse)
library(dplyr)
library(tidyr)
library(plyr)
library(readr)
library(readxl)
library(tibble)
library(ggplot2)
library(dada2); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library(phangorn)
library(DECIPHER)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(DESeq2)

# Set seed
set.seed(666)
setwd("C:/Users/ava135/OneDrive - University of Canterbury/Projects/Passenger_risk_incursions/r-data-analysis/16S")

asv_mat <- read_tsv("ASVs_counts_16S.tsv")
tax_mat <- read.csv("merged_tax.csv", header = TRUE)
samples_df <- read_excel("metadata.xlsx")

colnames(asv_mat)[1] <- "ASV"
colnames(tax_mat)[1] <- "ASV"

# Set rownames
asv_mat <- asv_mat %>% column_to_rownames("ASV")
tax_mat <- tax_mat %>% column_to_rownames("ASV")
samples_df <- samples_df %>% column_to_rownames("sample")

# Convert to matrices
asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

# Build phyloseq
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(samples_df)
)

ps   # 7197 taxa, 65 samples

######## filtering data step #################################################################################
##############################################################################################################

##remove uncgaracterised taxa
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) 
ps0 #7012 taxa and 65 samples
ps0 <- ps 

#Remove plant contaminations
ps_clean <- subset_taxa(ps, Phylum != "Streptophyta")
ps_clean #7012 taxa; 65 samples 

#Remove archaea contaminations
ps_clean <- subset_taxa(ps_clean, Phylum != "Methanobacteriota")
ps_clean #7007
ps <- ps_clean

######## decontam #################################################################################
##############################################################################################################

sample_names(ps) ##check to make sure the samples are all correct
lookup <- samples_df$`sample-name`
names(lookup) <- rownames(samples_df)
sample_names(ps) <- lookup[sample_names(ps)] ##fixed the sample names so they line up to the proper samples and blanks
sample_names(ps) ##chek


library(decontam)

sample_data(ps)$is.neg <- sample_names(ps) %in% c(
  "AS_T0_16S",
  "AS_T1_16S",
  "AS_T2_16S",
  "swabBLANK2_I6S",
  "SB2_16S",
  "16S_BLANK",
  "16S_BLANK_2_16S",
  "airBLANK_16S",
  "extBLANK_16S",
  "extBLANK_2_16S",
  "PCR1-B1-16S",
  "PCR1-B2-16S",
  "swabBLANK1_I6S",
  "XB_16S-1",
  "XB_16S-2",
  "AS_B1_16S",
  "XB_16S-4",
  "XB_16S-5",
  "PCR1-B1-I6S-1",
  "PCR1-B1-I6S-2",
  "SB1_16S"
)
                                                 
contamdf <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf$contaminant) #false = 6929; true = 78
contaminants <- rownames(contamdf)[contamdf$contaminant == TRUE]

ps_clean <- prune_taxa(!taxa_names(ps) %in% contaminants, ps)
ps_clean #6929 taxa; 65 samples

#remove taxa with 0 reads  and samples with <500 reads
ps_nozero <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
ps_filtered <- prune_samples(sample_sums(ps_nozero) >= 500, ps_nozero)
ps <- ps_filtered
ps <- ps_nozero
ps #6929


######## remove controls #################################################################################
##############################################################################################################

neg_samples <- c(
  "AS_T0_16S",
  "AS_T1_16S",
  "swabBLANK2_I6S",
  "SB2_16S",
  "16S_BLANK",
  "16S_BLANK_2_16S",
  "airBLANK_16S",
  "extBLANK_16S",
  "extBLANK_2_16S",
  "PCR1-B1-16S",
  "PCR1-B2-16S",
  "swabBLANK1_I6S",
  "XB_16S-1",
  "XB_16S-2",
  "AS_B1_16S",
  "XB_16S-4",
  "XB_16S-5",
  "PCR1-B1-I6S-1",
  "PCR1-B1-I6S-2",
  "SB1_16S",
  "AS_T2_16S"
)

ps_final <- subset_samples(ps, !sample_names(ps) %in% neg_samples)

ps_final #6929 taxa; 53 samples
sample_names(ps_final)

######## normalise data #################################################################################
##############################################################################################################

library(MicrobiomeStat)

otu_mat <- as(otu_table(ps_final), "matrix")
samples_df <- as(sample_data(ps_final), "data.frame")
tax_mat <- as(tax_table(ps_final), "matrix")
otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")

####################################################################################################
##################################DAT ANALYSIS SECTION##############################################
####################################################################################################

##used this for top 20 and top 15
top20 <- names(sort(taxa_sums(ps_final), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_final, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="sample.name", fill="Genus") + facet_wrap(~type, scales="free_x")


plot_bar(ps.top20, x = "sample.name", fill = "Genus") +
  facet_wrap(~type, scales = "free_x") +  # one panel per type
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Tent Number / Sample", y = "Relative Abundance", fill = "Genus")

##because top 20 are mutiple in same genus i plotted by Asv also

plot_bar(ps.top20, x = "sample.name", fill = "OTU") +  # or fill = "taxa_names"
  facet_wrap(~type, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Tent Number / Sample", y = "Relative Abundance", fill = "ASV")

#######sample stats and reads##############################################################
#######################################################################################

# Extract sample sums
sample_stats <- data.frame(
  sample = sample_names(ps_final),
  total_reads = sample_sums(ps_final)
)

# Add metadata
metadata <- data.frame(sample_data(ps_final))
metadata <- tibble::rownames_to_column(data.frame(sample_data(ps_final)), var = "sample")
sample_stats <- left_join(sample_stats, metadata, by = "sample")

# View the table
head(sample_stats)
write.csv(sample_stats, file = "sample_stats.csv", row.names = FALSE)


#######bray curtis ordination##############################################################
#######################################################################################

ps.clean <- prune_samples(sample_sums(ps_final) > 0, ps_final)
ps.rel <- transform_sample_counts(ps.clean, function(x) x / sum(x))
taxa_sums(ps.rel)
ps.rel <- prune_taxa(taxa_sums(ps.rel) > 0, ps.rel)

ord <- ordinate(ps.rel, method="PCoA", distance="bray")
library(ggrepel)

p <- plot_ordination(ps.rel, ord, color = "location...2", shape = "location...7") +
  geom_point(size = 4) +
  geom_text(aes(label = sample.name), vjust = -1, size = 3) +
  theme_classic()
p


#######alpha diversity##############################################################
#######################################################################################
# Remove zero-read samples
ps.clean <- prune_samples(sample_sums(ps_final) > 0, ps_final)
sample_sums(ps.clean)

# Find minimum sample depth
min_reads <- min(sample_sums(ps.clean))

# Rarefy to that depth
ps.rarefy <- prune_samples(sample_sums(ps.clean) >= 1000, ps.clean)

min_reads <- min(sample_sums(ps.rarefy))  # now > 1000
set.seed(123)
ps.rarefied <- rarefy_even_depth(ps.rarefy,
                                 sample.size = min_reads,
                                 rngseed = 123,
                                 verbose = FALSE)

# Calculate alpha diversity
alpha_div <- estimate_richness(ps.rarefied,
                               measures = c("Observed", "Shannon", "Simpson"))

# Add sample names as a column
alpha_div$sample <- rownames(alpha_div)
# Remove X from alpha_div sample names
alpha_div$sample <- sub("^X", "", alpha_div$sample)

# Then merge
alpha_div <- merge(alpha_div, meta_df,
                   by.x = "sample", by.y = "sample.name",
                   all.x = TRUE)

# Check
head(alpha_div[, c("sample", "Shannon", "type", "is.neg", "tent.number")])
# Extract metadata and merge correctly
meta_df <- data.frame(sample_data(ps.rarefied))
alpha_div <- merge(alpha_div, meta_df,
                   by.x = "sample", by.y = "sample.name",
                   all.x = TRUE)

# Check
head(alpha_div)
# Remove rows with NA in location...7.x
alpha_div_plot <- alpha_div %>% filter(!is.na(location...7.x))

##nz versus overseas
ggplot(alpha_div, aes(x = location...7.x, y = Shannon, fill = location...7.x)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic() +
  labs(x = "Overseas", y = "Shannon Diversity")

##nz versus overseasonly from tent
# Filter for tent samples only
alpha_div_tent <- alpha_div %>% 
  filter(type.x == "tent" & !is.na(location...7.x))  # also remove NAs in location

# Plot NZ vs Overseas for tent samples
ggplot(alpha_div_tent, aes(x = location...7.x, y = Shannon, fill = location...7.x)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic() +
  labs(x = "Overseas", y = "Shannon Diversity")

##by sample type
ggplot(alpha_div, aes(x = type.x, y = Shannon, fill = type.x)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic() +
  labs(x = "Sample Type", y = "Shannon Diversity")
