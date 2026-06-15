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
setwd("C:/Users/VaughanA/OneDrive - MWLR/Projects/microbial_hitchikers/R-analysis-mar2025/ITS")

asv_mat <- read_tsv("ASVs_counts_ITS.tsv")
tax_mat_ITS <- read.csv("merged_tax.csv", header = TRUE)
samples_df <- read_excel("metadata.xlsx")

colnames(asv_mat)[1] <- "ASV"
colnames(tax_mat_ITS)[1] <- "ASV"

# Set rownames
asv_mat <- asv_mat %>% column_to_rownames("ASV")
tax_mat_ITS <- tax_mat_ITS %>% column_to_rownames("ASV")
samples_df <- samples_df %>% column_to_rownames("seq_name")

# Convert to matrices
asv_mat <- as.matrix(asv_mat)
tax_mat_ITS <- as.matrix(tax_mat_ITS)

##check all samples are there and accounted for for sanity
asv_samples <- colnames(asv_mat)
meta_samples <- rownames(samples_df)
identical(asv_samples, meta_samples)

# In ASV matrix but not metadata
setdiff(asv_samples, meta_samples) ## 0
# In metadata but not ASV matrix
setdiff(meta_samples, asv_samples) ## 13 samples, those are the ones which have 0 reads

# Build phyloseq
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat_ITS),
  sample_data(samples_df)
)
ps   # 3358; 63



######## filtering data step #################################################################################
##############################################################################################################

##remove uncgaracterised taxa
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) 
ps0 #3244 taxa and 63 samples
ps0 <- ps 

##filter human contaminiation
ps_nohomo <- subset_taxa(ps, Family != "Hominidae")
ps_nohomo #2819, 62 samples

#Remove plant contaminations
ps_clean <- subset_taxa(ps_nohomo, Phylum != "Streptophyta")
ps_clean #2819 taxa; 63 samples 
ps <- ps_clean


######## decontam #################################################################################
##############################################################################################################

sample_names(ps) ##check to make sure the samples are all correct
lookup <- samples_df$`sample_name`
names(lookup) <- rownames(samples_df)
sample_names(ps) <- lookup[sample_names(ps)] ##fixed the sample names so they line up to the proper samples and blanks
sample_names(ps) ##chek


library(decontam)

sample_data(ps)$is.neg <- sample_names(ps) %in% c("ITS_BLANK",
                                                  "airBLANK_ITS",
                                                  "extBLANK_ITS",
                                                  "AS_B1_ITS",
                                                  "AS_T1_ITS",
                                                  "AS_T2_ITS",
                                                  "XB_ITS-4",
                                                  "XB_ITS-1",
                                                  "ITS_BLANK",
                                                  "ITS_BLANK_2",
                                                  "ITS_BLANK_3",
                                                  "PCR1-B1-ITS",
                                                  "PCR1-B2-ITS",
                                                  "swabBLANK1_ITS",
                                                  "swabBLANK2_ITS",
                                                  "XB_ITS-2",
                                                  "XB_ITS-3",
                                                  "PCR1-B1-I6S-1",
                                                  "PCR1-B1-I6S-2",
                                                  "PCR1-B1-ITS",
                                                  "PCR1-B2-ITS",
                                                  "SB1_ITS",
                                                  "SB2_ITS"
)

contamdf <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf$contaminant) #false = 2798; true = 21
contaminants <- rownames(contamdf)[contamdf$contaminant == TRUE]

ps_clean <- prune_taxa(!taxa_names(ps) %in% contaminants, ps)
ps_clean #2798 taxa; 62 samples

#remove taxa with 0 reads 
#ps_nozero <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
#ps <- ps_nozero
#ps #6929


######## remove controls #################################################################################
##############################################################################################################

neg_samples <- c("ITS_BLANK",
                 "airBLANK_ITS",
                 "extBLANK_ITS",
                 "AS_B1_ITS",
                 "AS_T1_ITS",
                 "AS_T2_ITS",
                 "XB_ITS-4",
                 "XB_ITS-1",
                 "ITS_BLANK_2",
                 "ITS_BLANK_3",
                 "PCR1-B1-ITS",
                 "PCR1-B2-ITS",
                 "swabBLANK1_ITS",
                 "swabBLANK2_ITS",
                 "XB_ITS-2",
                 "XB_ITS-3",
                 "PCR1-B1-I6S-1",
                 "PCR1-B1-I6S-2",
                 "PCR1-B1-ITS",
                 "PCR1-B2-ITS",
                 "SB1_ITS",
                 "SB2_ITS"
)

ps_final <- subset_samples(ps, !sample_names(ps) %in% neg_samples)

ps_final #2819 taxa; 52 samples
sample_names(ps_final)

##remove samples that have extremly low reads
ps_final <- prune_samples(sample_sums(ps_final) > 100, ps_final)
sample_sums(ps_final)
ps_final #2819; 39 samples
ps_final_ITS <- ps_final
#############################################################
################### Normalize data ##########################
#############################################################

##check format of phyloseq object is correct
phyloseq::taxa_are_rows(ps_final_ITS)

data.obj.ITS <- mStat_convert_phyloseq_to_data_obj(ps_final_ITS)
norm.data.obj <- mStat_normalize_data(data.obj.ITS, "TSS")$data.obj.norm
mStat_validate_data(norm.data.obj) ##validation pass

mStat_rarefy_data(data.obj=ps_final, depth = 500)


####################################################################################################
##################################DAT ANALYSIS SECTION##############################################
####################################################################################################


##see extra R script for analysis




