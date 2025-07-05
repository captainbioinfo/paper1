# Load required packages
library(tidyverse)
library(RColorBrewer)   # for 'Dark2' palette
df <- read_csv("phylon_only.csv")  # Replace with your actual file path
# Reshape to long format
df_long <- df %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Count")
ggplot(df_long, aes(x = Sample, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Genus Abundance Across Samples",
       x = "Sample", y = "Read Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Calculate relative abundance
df_relative <- df_long %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Count / sum(Count) * 100)

# Create the plot with a new colour scheme
p <- ggplot(df_relative, aes(x = Sample, y = RelAbundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  theme_classic(base_size = 22) +            # white background
  scale_fill_brewer(palette = "Dark2") +     # << new colour palette
  labs(
    title = "  ",
    x = "Studies (Bovine Mastitis)",
    y = "Relative Abundance (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title  = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 24, face = "bold", hjust = 0.5)  # centred
  )
plot(p)
# Save at high resolution with white background
ggsave("Relative_Abundance_Barplot23.png",  plot = p,
       width = 10, height = 12, dpi = 600, bg = "white")

# Optional TIFF for journal submission
ggsave("Relative_Abundance_Barplot23.tiff", plot = p,
       width = 10, height = 12, dpi = 600, bg = "white",
       compression = "lzw")
