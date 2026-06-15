#############################
## RELATIVE ABUNDANCE PLOT 
############################

# Filter out unknown/ambiguous classes
ps_class <- tax_glom(ps_final_ITS, taxrank = "Class")
ps_class <- subset_taxa(ps_class, !grepl("Incertae_sedis|unassigned|Unknown", Class))
ps_class
# Step 2: Keep top 15 most abundant classes
top15 <- names(sort(taxa_sums(ps_class), decreasing = TRUE))[1:15]
ps_top15_class <- prune_taxa(top15, ps_class)

#make relative 
ps_top15_class <- transform_sample_counts(ps_top15_class, function(x) x / sum(x))

# Step 3: Melt to dataframe
df <- psmelt(ps_top15_class)

#Sort out the colours

myPal <- tax_palette(
  data = ps_class, rank = "Class", n = 15, pal = "kelly",
  add = c(Other = "white")
)
tax_palette_plot(myPal) ##brewerPlus 16S and kelly ITS

# Common theme
base_theme <- theme_minimal(base_size = 11, base_family = "Montserrat") +
  theme(
    text = element_text(color = "black"),  # <-- this line makes all text black
    axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.margin = margin(5, -5, 5, -5),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  )
# Step 5: Subset by location
df_overseas  <- subset(df, overseas_prev_6months == "Yes")
df_nz  <- subset(df, overseas_prev_6months == "no")


# Remove y-axis text for all but first plot
no_y_axis <- theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank()
)

# Plot 1 — keep y-axis
p1 <- ggplot(df_overseas, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "Overseas", y = "Relative abundance", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + theme(plot.title = element_text(face = "bold.italic"))

# Plot 2
p2 <- ggplot(df_nz, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "New Zealand", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + no_y_axis + theme(plot.title = element_text(face = "bold.italic"))

# Combine with adjusted widths & shared legend
final_plot <- (p2+p1) +
  plot_layout(guides = "collect", widths = c(2, 1)) & 
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", family = "Montserrat", size = 13),
    legend.text = element_text(size = 12),
    legend.box = "vertical",
    legend.title.align = 0.5
  ) &
  guides(fill = guide_legend(title = "Class", title.position = "top"))

# Show plot
final_plot


# Aggregate to Genus level
ps_genus <- tax_glom(ps_final, taxrank = "Genus")

# Remove ambiguous genera
ps_genus <- subset_taxa(
  ps_genus,
  !grepl("Incertae_sedis|unassigned|Unknown|uncultured", Genus)
)

# Keep top 50 most abundant genera
top50 <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:50]

ps_top50_genus <- prune_taxa(top50, ps_genus)

# Convert to relative abundance
ps_top50_genus <- transform_sample_counts(
  ps_top50_genus,
  function(x) x / sum(x)
)

# Convert to dataframe
df_genus <- psmelt(ps_top50_genus)

df_overseas <- subset(
  df_genus,
  overseas_prev_6months == "Yes"
)

df_nz <- subset(
  df_genus,
  overseas_prev_6months == "no"
)

# Generate 50 colours
cols50 <- colorRampPalette(
  brewer.pal(12, "Paired")
)(50)

base_theme <- theme_minimal(
  base_size = 11,
  base_family = "Montserrat"
) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 90,
      size = 5,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.margin = margin(5, -5, 5, -5),
    plot.title = element_text(
      face = "bold.italic",
      size = 13,
      hjust = 0.5
    )
  )

# Remove y-axis for second plot
no_y_axis <- theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank()
)

# Overseas plot
p1 <- ggplot(
  df_overseas,
  aes(x = Sample, y = Abundance, fill = Genus)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols50) +
  labs(
    title = "Overseas",
    y = "Relative abundance",
    x = NULL
  ) +
  base_theme

# New Zealand plot
p2 <- ggplot(
  df_nz,
  aes(x = Sample, y = Abundance, fill = Genus)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols50) +
  labs(
    title = "New Zealand",
    x = NULL
  ) +
  base_theme +
  no_y_axis

################################################################################
## COMBINE
################################################################################

final_plot <- (p2 + p1) +
  plot_layout(
    guides = "collect",
    widths = c(2, 1)
  ) &
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      family = "Montserrat",
      size = 13
    ),
    legend.text = element_text(size = 8),
    legend.box = "vertical",
    legend.title.align = 0.5
  ) &
  guides(
    fill = guide_legend(
      title = "Genus",
      title.position = "top",
      ncol = 2
    )
  )

# Show plot
final_plot


###compare overseas and nz for genera presence/ absence

# Genera present in NZ
nz_genera <- unique(
  df_genus$Genus[
    df_genus$overseas_prev_6months == "no" &
      df_genus$Abundance > 0
  ]
)

# Genera present in Overseas
overseas_genera <- unique(
  df_genus$Genus[
    df_genus$overseas_prev_6months == "Yes" &
      df_genus$Abundance > 0
  ]
)

# NZ only
nz_only <- setdiff(nz_genera, overseas_genera)

# Overseas only
overseas_only <- setdiff(overseas_genera, nz_genera)

nz_only
overseas_only
