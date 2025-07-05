# Install (if not already installed)
install.packages("RColorBrewer")
# Install if not already installed
install.packages("ggrepel")

# Load the library
library(ggrepel)

# Load the package
library(RColorBrewer)

# Load necessary libraries
library(vegan)     # for distance and ordination
library(ggplot2)   # for plotting
library(ape)       # for PCoA (pcoa function)
library(reshape2)  # for data manipulation (optional)

# Load your data
data <- read.csv("PRJEN43443_edited.csv", row.names = 1, check.names = FALSE)

data_t <- t(data)
data_t <- as.data.frame(data_t)
data_t <- data_t / rowSums(data_t)  # Normalize

# Bray-Curtis distance
bc_dist <- vegdist(data_t, method = "bray")

# PCoA
pcoa_res <- pcoa(bc_dist)
pcoa_df <- as.data.frame(pcoa_res$vectors)
pcoa_df$Sample <- rownames(pcoa_df)

# Optional: Add sample groups manually or load metadata
# Example (replace with real metadata):
# pcoa_df$Group <- factor(c(rep("Group1", 10), rep("Group2", 10), rep("Group3", 10)))

# Colors (you can change 'Group' to actual metadata column if available)
set.seed(42)
n_colors <- length(unique(pcoa_df$Sample))  # Or pcoa_df$Group if available
colors <- colorRampPalette(brewer.pal(8, "Set1"))(n_colors)

# Plot
p <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, label = Sample)) +
  geom_point(aes(color = Sample), size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20, segment.size = 0.2) +
  scale_color_manual(values = colors) +
  xlab(paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)")) +
  ylab(paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_minimal(base_size = 16) +
  ggtitle("Improved PCoA Plot (Bray-Curtis Distance)") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.8)
  )
plot(p)

ggsave("PCoA_plot.png", plot = p, width = 8, height = 6, dpi = 300)
ggsave("PCoA_plot_white_bg.png", plot = p, width = 8, height = 6, dpi = 300, bg = "white")
# 
