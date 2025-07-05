# Install (if not already installed)
install.packages("RColorBrewer")
install.packages("ggrepel")

# Load libraries
library(ggrepel)
library(RColorBrewer)
library(ggplot2)

# Load your data
data <- read.csv("PRJEN43443_edited.csv", row.names = 1, check.names = FALSE)

# Transpose and normalize
data_t <- t(data)
data_t <- as.data.frame(data_t)
data_t <- data_t / rowSums(data_t)  # Normalize

# Perform PCA
pca_res <- prcomp(data_t, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)  # Extract PCA components
pca_df$Sample <- rownames(pca_df)

# Set colors
set.seed(42)
n_colors <- length(unique(pca_df$Sample))
colors <- colorRampPalette(brewer.pal(8, "Set1"))(n_colors)

# Plot PCA
p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(aes(color = Sample), size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20, segment.size = 0.2) +
  scale_color_manual(values = colors) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal(base_size = 16) +
  ggtitle("Improved PCA Plot") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.8)
  )

plot(p)

# Save the plot with a white background
ggsave("PCA_plot_white_bg.png", plot = p, width = 8, height = 6, dpi = 300, bg = "white")
