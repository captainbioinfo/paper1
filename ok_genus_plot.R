library(tidyverse)
library(RColorBrewer)

# Read data
df_genus <- read_csv("genus_only.csv")

# Reshape to long format
df_genus_long <- df_genus %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Count")

# Calculate relative abundance
df_genus_rel <- df_genus_long %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Count / sum(Count) * 100)

# Generate a long color palette
n_colors <- length(unique(df_genus_rel$Genus))
colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)

# Create the plot
p_genus <- ggplot(df_genus_rel, aes(x = Sample, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = colors) +
  labs(
    title = " ",
    x = "Studies (Bovine Mastitis)", y = "Relative Abundance (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.25)
  )
plot(p_genus)
# Save as high-resolution image
ggsave("Genus_Relative_Abundance22.png", plot = p_genus,
       width = 20, height = 14, dpi = 600, bg = "white")

ggsave("Genus_Relative_Abundance22.tiff", plot = p_genus,
       width = 20, height = 14, dpi = 600, bg = "white", compression = "lzw")
##################################################################################################################

