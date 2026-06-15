#############################################################
####################### Clamtest ## #########################
#############################################################

library(phyloseq)
library(vegan)
library(labdsv)

##subset into travel history
overseas <- subset_samples(ps_final_16S, overseas_prev_6months == "Yes")
nz       <- subset_samples(ps_final_16S, overseas_prev_6months == "no")


otu_overseas <- as(otu_table(overseas), "matrix")
otu_nz       <- as(otu_table(nz), "matrix")

# Fix orientation (samples must be rows)
if (taxa_are_rows(overseas)) otu_overseas <- t(otu_overseas)
if (taxa_are_rows(nz)) otu_nz <- t(otu_nz)


otu_all <- rbind(otu_overseas, otu_nz)

group <- c(
  rep("Overseas", nrow(otu_overseas)),
  rep("New Zealand", nrow(otu_nz))
)

##clean the data
# remove taxa with zero total abundance
otu_all <- otu_all[, colSums(otu_all) > 0]

##CLAM test
clam_res <- clamtest(otu_all, group)


cat("Total taxa:", length(clam_res), "\n")


clam_df <- as.data.frame(clam_res)

nz_specialists <- subset(clam_df, Classes == "Specialist_New Zealand")
overseas_specialists <- subset(clam_df, Classes == "Specialist_Overseas")
generalists <- subset(clam_df, Classes == "Generalist")

write.csv(nz_specialists, "CLAM_NZ_specialists.csv", row.names = FALSE)
write.csv(overseas_specialists, "CLAM_Overseas_specialists.csv", row.names = FALSE)
write.csv(generalists, "CLAM_Generalists.csv", row.names = FALSE)

write.csv(clam_df, "CLAM_all_results.csv", row.names = FALSE)

table(clam_df$Classes)
prop.table(table(clam_df$Classes))
clam_summary <- as.data.frame(table(clam_df$Classes))
colnames(clam_summary) <- c("Category", "Count")
clam_summary$Proportion <- clam_summary$Count / sum(clam_summary$Count)

clam_summary
write.csv(clam_summary, "CLAM_summary_ASV_level.csv", row.names = FALSE)

ps_class <- tax_glom(ps_final, taxrank = "Class")
otu_class <- as(otu_table(ps_class), "matrix")

if (taxa_are_rows(ps_class)) {
  otu_class <- t(otu_class)
}
overseas_class <- subset_samples(ps_class, overseas_prev_6months == "Yes")
nz_class <- subset_samples(ps_class, overseas_prev_6months == "no")

otu_overseas <- as(otu_table(overseas_class), "matrix")
otu_nz <- as(otu_table(nz_class), "matrix")

if (taxa_are_rows(overseas_class)) otu_overseas <- t(otu_overseas)
if (taxa_are_rows(nz_class)) otu_nz <- t(otu_nz)

otu_all_class <- rbind(otu_overseas, otu_nz)

group_class <- c(
  rep("Overseas", nrow(otu_overseas)),
  rep("New Zealand", nrow(otu_nz))
)

otu_all_class <- otu_all_class[, colSums(otu_all_class) > 0]
stopifnot(nrow(otu_all_class) == length(group_class))

clam_class <- clamtest(otu_all_class, group_class)
clam_class_df <- as.data.frame(clam_class)

class_summary <- as.data.frame(table(clam_class_df$Classes))
colnames(class_summary) <- c("Category", "Count")
class_summary$Proportion <- class_summary$Count / sum(class_summary$Count)

class_summary
write.csv(clam_class_df, "CLAM_class_level_results.csv", row.names = FALSE)
write.csv(class_summary, "CLAM_class_level_summary.csv", row.names = FALSE)

########################
nz_specialists
nz_taxonomy <- tax_df[tax_df$ASV %in% nz_specialists, ]
nz_taxonomy
nz_taxonomy[, c("ASV", colnames(nz_taxonomy)[colnames(nz_taxonomy) != "ASV"])]

nz_ids <- trimws(as.character(nz_specialists$Species))
tax_df <- as.data.frame(tax_mat)
tax_df$ASV <- trimws(rownames(tax_df))
nz_taxonomy <- tax_df[tax_df$ASV %in% nz_ids, ]

################# do for genus level
ps_genus <- tax_glom(ps_final, taxrank = "Genus")

# remove empty/unknown genera
ps_genus <- subset_taxa(
  ps_genus,
  !is.na(Genus) &
    Genus != "" &
    !grepl("unclassified|uncultured|Incertae_sedis|Unknown", Genus)
)


otu_genus <- as(otu_table(ps_genus), "matrix")

if (taxa_are_rows(ps_genus)) {
  otu_genus <- t(otu_genus)
}


group <- sample_data(ps_genus)$overseas_prev_6months
group <- ifelse(group == "Yes", "Overseas", "New Zealand")

##clean dataset
otu_genus <- otu_genus[, colSums(otu_genus) > 0]

stopifnot(nrow(otu_genus) == length(group))


clam_genus <- clamtest(otu_genus, group)
clam_genus_df <- as.data.frame(clam_genus)


genus_summary <- as.data.frame(table(clam_genus_df$Classes))
colnames(genus_summary) <- c("Category", "Count")
genus_summary$Proportion <- genus_summary$Count / sum(genus_summary$Count)

print(genus_summary)


write.csv(clam_genus_df, "CLAM_genus_all_results.csv", row.names = FALSE)
write.csv(genus_summary, "CLAM_genus_summary.csv", row.names = FALSE)

write.csv(subset(clam_genus_df, Classes == "Generalist"),
          "CLAM_genus_generalists.csv", row.names = FALSE)

write.csv(subset(clam_genus_df, Classes == "Specialist_New Zealand"),
          "CLAM_genus_NZ_specialists.csv", row.names = FALSE)

write.csv(subset(clam_genus_df, Classes == "Specialist_Overseas"),
          "CLAM_genus_overseas_specialists.csv", row.names = FALSE)

###double check the resykts
table(clam_genus_df$Classes)

make_clam_summary <- function(clam_obj, level_name) {
  df <- as.data.frame(clam_obj)
  
  summary <- as.data.frame(table(df$Classes))
  colnames(summary) <- c("Category", "Count")
  summary$Proportion <- summary$Count / sum(summary$Count)
  
  summary$Level <- level_name
  return(summary)
}

asv_summary <- make_clam_summary(clam_res, "ASV")
genus_summary <- make_clam_summary(clam_genus, "Genus")
class_summary <- make_clam_summary(clam_genus, "Class")

clam_summary_all <- rbind(
  asv_summary,
  genus_summary,
  class_summary
)

clam_summary_all

get_clam_taxa <- function(clam_obj, level_name) {
  df <- as.data.frame(clam_obj)
  
  taxa_list <- split(df$Species, df$Classes)
  
  out <- data.frame(
    Category = names(taxa_list),
    Taxa = sapply(taxa_list, function(x) paste(x, collapse = ", ")),
    Level = level_name
  )
  
  rownames(out) <- NULL
  return(out)
}

asv_taxa <- get_clam_taxa(clam_res, "ASV")
genus_taxa <- get_clam_taxa(clam_genus, "Genus")
class_taxa <- get_clam_taxa(clam_class, "Class")
