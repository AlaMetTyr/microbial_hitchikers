library(ComplexUpset)
library(ggplot2)

###more robust clam with proper definition  of biome membership using prevalence threshold


ps_genus <- tax_glom(ps_final_ITS, taxrank = "Genus")

ps_genus <- subset_taxa(
  ps_genus,
  !is.na(Genus) &
    Genus != "" &
    !grepl("uncultured|unclassified|Unknown|Incertae_sedis", Genus)
)

##convert to matrix
otu <- as(otu_table(ps_genus), "matrix")

if (taxa_are_rows(ps_genus)) {
  otu <- t(otu)
}

meta <- as.data.frame(sample_data(ps_genus))

otu_pa <- otu
otu_pa[otu_pa > 0] <- 1

###define biome membership
biomes <- unique(meta$biome)

prev_threshold <- 0.05  # adjust: 0.05–0.2 depending on strictness

genus_lists <- lapply(biomes, function(b) {
  
  samples_b <- rownames(meta[meta$biome == b, , drop = FALSE])
  
  if (length(samples_b) == 0) return(character(0))
  
  prev <- colMeans(otu_pa[samples_b, , drop = FALSE])
  
  names(prev[prev >= prev_threshold])
})

names(genus_lists) <- biomes

###
all_taxa <- unique(unlist(genus_lists))

upset_df <- data.frame(Taxon = all_taxa)

for (b in biomes) {
  upset_df[[b]] <- all_taxa %in% genus_lists[[b]]
}

###Plot
ComplexUpset::upset(
  upset_df,
  intersect = biomes,
  name = "Biomes",
  width_ratio = 0.2
)

