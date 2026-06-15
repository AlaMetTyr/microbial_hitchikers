########relationship between asv and last use

# Pull sample metadata
meta <- data.frame(sample_data(ps))
meta$sample <- rownames(meta)  # ensure sample names are a column
alpha_div$sample <- rownames(alpha_div)
library(dplyr)

merged_df <- alpha_div %>%
  left_join(meta, by="sample")
library(dplyr)
head(merged_df)

# Convert character to numeric
merged_df$days_between.num <- as.numeric(merged_df$days_between.x)

# Check conversion
str(merged_df$days_between.num)
summary(merged_df$days_between.num)

# Pearson correlation
cor.test(merged_df$Observed, merged_df$days_between.num)

# Or Spearman if the data is not normal
cor.test(merged_df$Observed, merged_df$days_between.num, method="spearman")


library(ggplot2)

ggplot(merged_df, aes(x=days_between.num, y=Observed)) +
  geom_point(size=3, alpha=0.7) +
  geom_smooth(method="lm", col="blue") +
  labs(x="Days since last use", y="Number of ASVs") +
  theme_minimal()
