# Sample coverage
#Create a sample sum table to look at coverage metrics
# Extract read counts per sample
sample_sum_df <- data.frame(
  SampleID = names(sample_sums(ps_final)),
  Reads = sample_sums(ps_final)
)

# Extract metadata (e.g., species) from phyloseq object
metadata <- data.frame(sample_data(ps_final))

# Combine read counts with species info
sample_sum_df <- cbind(sample_sum_df, metadata[rownames(metadata) %in% sample_sum_df$SampleID, ])


# Plot histogram
ggplot(sample_sum_df, aes(x = Reads)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Save to CSV
write.csv(sample_sum_df, "ITS_Reads_Per_Sample.csv", row.names = FALSE)
