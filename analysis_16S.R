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
library(decontam)
library(MicrobiomeStat)
library(microbiome)
library(microViz)
library(ggplot2)
library(patchwork)
library(RColorBrewer)


# Set seed
set.seed(666)
setwd("C:/Users/VaughanA/OneDrive - MWLR/Projects/microbial_hitchikers/R-analysis-mar2025/16S")

asv_mat <- read_tsv("ASVs_counts_16S.tsv")
tax_mat_16S <- read.csv("merged_tax.csv", header = TRUE)
samples_df <- read_excel("metadata.xlsx")

colnames(asv_mat)[1] <- "ASV"
colnames(tax_mat_16S)[1] <- "ASV"

# Set rownames
asv_mat <- asv_mat %>% column_to_rownames("ASV")
tax_mat_16S <- tax_mat_16S %>% column_to_rownames("ASV")
samples_df <- samples_df %>% column_to_rownames("seq_name")

# Convert to matrices
asv_mat <- as.matrix(asv_mat)
tax_mat_16S <- as.matrix(tax_mat_16S)

##check all samples are there and accounted for for sanity
asv_samples <- colnames(asv_mat)
meta_samples <- rownames(samples_df)
identical(asv_samples, meta_samples)

# In ASV matrix but not metadata
setdiff(asv_samples, meta_samples) ## 0
# In metadata but not ASV matrix
setdiff(meta_samples, asv_samples) ## 11 samples, those are the ones which have 0 reads

# Build phyloseq
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat_16S),
  sample_data(samples_df)
)

ps   # 7197 taxa, 65 samples

#############################################################
######## filtering data step ################################
#############################################################

##remove uncgaracterised taxa
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) 
ps0 #7012 taxa and 65 samples
ps0 <- ps 

##filter human contaminiation
ps_nohomo <- subset_taxa(ps, Family != "Hominidae")
ps_nohomo #5735, 65 samples

#Remove plant contaminations
ps_clean <- subset_taxa(ps_nohomo, Phylum != "Streptophyta")
ps_clean #5735 taxa; 65 samples 

#Remove archaea contaminations
ps_clean <- subset_taxa(ps_clean, Phylum != "Methanobacteriota")
ps_clean #5730; 65 samples
ps <- ps_clean

#############################################################
################### DECONTAM ################################
#############################################################

sample_names(ps) ##check to make sure the samples are all correct
lookup <- samples_df$`sample_name`
names(lookup) <- rownames(samples_df)
sample_names(ps) <- lookup[sample_names(ps)] ##fixed the sample names so they line up to the proper samples and blanks
sample_names(ps) ##chek

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
table(contamdf$contaminant) #false = 5658; true = 72
contaminants <- rownames(contamdf)[contamdf$contaminant == TRUE]

ps_clean <- prune_taxa(!taxa_names(ps) %in% contaminants, ps)
ps_clean #5658 taxa; 65 samples

#remove taxa with 0 reads 
#ps_nozero <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
#ps <- ps_nozero
#ps #6929

#############################################################
################### Remove controls #########################
#############################################################

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

ps_final #5730 taxa; 53 samples
sample_names(ps_final)
ps_final_16S <- ps_final
##remove samplesps_final_16S##remove samples that have extremly low reads
ps_final_16S <- prune_samples(sample_sums(ps_final_16S) > 100, ps_final)
sample_sums(ps_final_16S)
ps_final_16S #5730; 38 samples
sample_sums(ps_final_16S)

#############################################################
################### Normalize data ##########################
#############################################################

##check format of phyloseq object is correct
phyloseq::taxa_are_rows(ps_final_16S)

data.obj.16S <- mStat_convert_phyloseq_to_data_obj(ps_final_16S)
norm.data.obj <- mStat_normalize_data(data.obj.16S, "TSS")$data.obj.norm
mStat_validate_data(norm.data.obj) ##validation pass

mStat_rarefy_data(data.obj=ps_final, depth = 500)


####################################################################################################
##################################DAT ANALYSIS SECTION##############################################
####################################################################################################


##see additional R script for full scripts








