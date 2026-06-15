#############################################################
################### Beta diversity #########################
#############################################################

library(aplot)
library(ggplotify)
library(dplyr)
library(ggplot2)
library(tibble)

#generate beta ordination based on bray-curtis dissimilarities
beta.ord <- generate_beta_ordination_single(
  data.obj   = norm.data.obj,
  dist.obj   = NULL,
  pc.obj     = NULL,
  subject.var = NULL,
  time.var   = NULL,
  t.level    = NULL,
  group.var  = "biome",
  strata.var = NULL,
  dist.name  = "BC",
  base.size  = 18,
  theme.choice = "bw",
  palette    = NULL,
  pdf        = FALSE
)
beta.ord


## bray curtis and pcoa

# Calculate Bray-Curtis distances
beta.bc <- mStat_calculate_beta_diversity(
  data.obj = norm.data.obj,
  dist.name = "BC"
)

bc.dist <- beta.bc$BC

# Run PCoA
pcoa_16s <- pcoa(bc.dist)

# Extract coordinates
pcoa_df <- as.data.frame(pcoa_16s$vectors) %>%
  rownames_to_column("Sample") %>%
  left_join(
    norm.data.obj$meta.dat %>% rownames_to_column("Sample"),
    by = "Sample"
  )

# Ensure grouping variable is clean factor
pcoa_df$overseas_prev_6months <- factor(
  pcoa_df$overseas_prev_6months,
  levels = c("no", "Yes")
)



var1 <- round(pcoa_16s$values$Relative_eig[1] * 100, 2)
var2 <- round(pcoa_16s$values$Relative_eig[2] * 100, 2)


overseas_palette <- c(
  "no" = "#2C7BB6",
  "Yes" = "#D7191C"
)



ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = overseas_prev_6months)) +
  geom_point(size = 6, alpha = 1) +
  stat_ellipse(size = 1, alpha = 1) +
  
  scale_color_manual(
    values = overseas_palette,
    name = "Overseas (last 6 months)",
    labels = c("no" = "No", "Yes" = "Yes")
  ) +
  
  coord_fixed(
    sqrt(pcoa_16s$values$Relative_eig[1] /
           pcoa_16s$values$Relative_eig[2])
  ) +
  
  labs(
    x = paste0("Axis 1 (", var1, "%)"),
    y = paste0("Axis 2 (", var2, "%)")
  ) +
  
  theme(
    legend.position = "right",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14, margin = margin(t = 10, r = 10)),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

##ggplot niceness for overseas
pcoa_df$overseas_prev_6months <- factor(
  pcoa_df$overseas_prev_6months,
  levels = c("no", "Yes"),
  labels = c("New Zealand", "Overseas")
)

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = overseas_prev_6months)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = overseas_prev_6months), linewidth = 1) +
  labs(
    x = paste0("Axis 1 (", var1, "%)"),
    y = paste0("Axis 2 (", var2, "%)"),
    colour = "Location of use"
  ) +
  theme_minimal()

##now do it to see if theres a trend by biome
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = biome)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = biome), linewidth = 1) +
  labs(
    x = paste0("Axis 1 (", var1, "%)"),
    y = paste0("Axis 2 (", var2, "%)"),
    colour = "Biome"
  ) +
  theme_minimal()


#PERMANOVA
library(GUniFrac)

permanova_res <- generate_beta_test_single(
  data.obj = norm.data.obj,
  dist.obj = NULL,
  time.var = NULL,
  t.level = NULL,
  group.var = "overseas_prev_6months",   # <-- CHANGE HERE
  adj.vars = "biome",  # optional covariate
  dist.name = c("BC")
)

otu <- norm.data.obj$feature.tab  # adjust if needed
otu <- t(otu)   # fix vegan oritentation
bc.dist <- vegdist(otu, method = "bray")
bd <- betadisper(bc.dist, norm.data.obj$meta.dat$overseas_prev_6months)
bd

permutest(bd, pairwise = TRUE)

bd_df <- data.frame(
  Biome = bd$group,
  DistanceToCentroid = bd$distances
)

ggplot(bd_df, aes(x = Biome, y = DistanceToCentroid)) +
  geom_boxplot(aes(fill = Biome), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(fill = Biome),
              color = "grey10",
              width = 0.2, size = 2.5,
              shape = 21, stroke = 0.8) +
  theme_minimal(base_family = "Montserrat") +
  labs(x = NULL, y = "Distance to centroid") +
  theme(legend.position = "none")

ps.clean <- prune_samples(sample_sums(norm.data.obj) > 0, norm.data.obj)
ps.rel <- transform_sample_counts(ps.clean, function(x) x / sum(x))

ps <- norm.data.obj$ps
ps.clean <- prune_samples(sample_sums(ps) > 0, ps)
ps.rel <- transform_sample_counts(ps.clean, function(x) x / sum(x))

ord <- ordinate(ps.rel, method = "PCoA", distance = "bray")

p <- plot_ordination(
  ps.rel,
  ord,
  color = "biome",
  shape = "overseas_prev_6months"
) +
  geom_point(size = 4) +
  theme_classic()

p
