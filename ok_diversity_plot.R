
# =============================================================================
# Microbial Diversity Analysis - Alpha and Beta Diversity Box Plots
# Publication-Ready R Script
# =============================================================================

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  vegan,          # For diversity calculations
  ggplot2,        # For plotting
  dplyr,          # For data manipulation
  tidyr,          # For data reshaping
  RColorBrewer,   # For color palettes
  gridExtra,      # For arranging plots
  cowplot,        # For publication-ready plots
  ggpubr,         # For statistical comparisons
  reshape2,       # For data melting
  scales,         # For plot scaling
  viridis,        # For color palettes
  ggrepel         # For text labels without overlap
)

# Set working directory (adjust as needed)
# setwd("your/working/directory")

# =============================================================================
# Data Input and Preparation
# =============================================================================

# Read CSV file
# Replace "your_data.csv" with the actual path to your CSV file
csv_file <- "nor_diversity_plot_data_cleaned.csv"  # Change this to your CSV file path

# Check if file exists
if (!file.exists(csv_file)) {
  stop(paste("CSV file not found:", csv_file, 
             "\nPlease ensure the file exists and update the csv_file variable with the correct path."))
}

# Read the CSV data
cat("Reading CSV file:", csv_file, "\n")
data_matrix <- read.csv(csv_file, 
                        header = TRUE, 
                        row.names = NULL,
                        stringsAsFactors = FALSE,
                        check.names = FALSE)  # Preserve original column names

# Display the first few rows to verify data loading
cat("First 5 rows of loaded data:\n")
print(head(data_matrix, 5))

# Check if first column contains taxa names
if (is.character(data_matrix[,1]) && !is.numeric(data_matrix[,1])) {
  # First column is taxa names
  taxa_names <- data_matrix[,1]
  abundance_data <- data_matrix[,-1]
} else {
  # Assume row names are taxa (if CSV has row names) or create generic names
  if (!is.null(rownames(data_matrix)) && rownames(data_matrix)[1] != "1") {
    taxa_names <- rownames(data_matrix)
    abundance_data <- data_matrix
  } else {
    # Create generic taxa names if none provided
    taxa_names <- paste("Taxa", 1:nrow(data_matrix), sep = "_")
    abundance_data <- data_matrix
  }
}

# Ensure all abundance data is numeric
abundance_data[] <- lapply(abundance_data, function(x) as.numeric(as.character(x)))

# Handle any NA values (convert to 0)
abundance_data[is.na(abundance_data)] <- 0

# Create clean data matrix
data_matrix <- data.frame(Taxa = taxa_names, abundance_data)

cat("\nData summary after processing:\n")
cat("Number of taxa:", nrow(data_matrix), "\n")
cat("Number of samples:", ncol(data_matrix) - 1, "\n")
cat("Sample names:", paste(names(data_matrix)[-1], collapse = ", "), "\n")

# Convert to matrix format (taxa as rows, samples as columns)
abundance_matrix <- as.matrix(data_matrix[, -1])
rownames(abundance_matrix) <- data_matrix$Taxa
colnames(abundance_matrix) <- names(data_matrix)[-1]

# Transpose for vegan package (samples as rows)
community_matrix <- t(abundance_matrix)

cat("\nFinal data dimensions:\n")
cat("Samples:", nrow(community_matrix), "\n")
cat("Taxa:", ncol(community_matrix), "\n")
cat("Sample names:", paste(rownames(community_matrix), collapse = ", "), "\n")

# =============================================================================
# Alpha Diversity Calculations
# =============================================================================

# Calculate multiple alpha diversity indices
alpha_diversity <- data.frame(
  Sample = rownames(community_matrix),
  Shannon = diversity(community_matrix, index = "shannon"),
  Simpson = diversity(community_matrix, index = "simpson"),
  InvSimpson = diversity(community_matrix, index = "invsimpson"),
  Richness = specnumber(community_matrix),
  Evenness = diversity(community_matrix, index = "shannon") / log(specnumber(community_matrix))
)

# Handle NaN values in Evenness (when richness = 0 or 1)
alpha_diversity$Evenness[is.nan(alpha_diversity$Evenness)] <- 0

print("Alpha Diversity Summary:")
print(summary(alpha_diversity[, -1]))

# =============================================================================
# Beta Diversity Calculations
# =============================================================================

# Calculate Bray-Curtis dissimilarity matrix
bray_curtis <- vegdist(community_matrix, method = "bray")

# Calculate Jaccard dissimilarity matrix
jaccard <- vegdist(community_matrix, method = "jaccard", binary = TRUE)

# Convert distance matrices to data frames for plotting
bray_df <- data.frame(
  Comparison = "Bray-Curtis",
  Dissimilarity = as.vector(bray_curtis)
)

jaccard_df <- data.frame(
  Comparison = "Jaccard",
  Dissimilarity = as.vector(jaccard)
)

# Combine beta diversity data
beta_diversity <- rbind(bray_df, jaccard_df)

print("Beta Diversity Summary:")
print(aggregate(Dissimilarity ~ Comparison, beta_diversity, summary))

# =============================================================================
# Data Preparation for Sample-wise Plotting
# =============================================================================

# Prepare alpha diversity data for sample-wise plotting
alpha_long <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson, InvSimpson, Evenness), 
               names_to = "Index", 
               values_to = "Value")

# Create beta diversity data for sample-wise analysis
# Convert distance matrix to long format with sample pairs
bray_matrix <- as.matrix(bray_curtis)
jaccard_matrix <- as.matrix(jaccard)

# Create sample-wise beta diversity (mean dissimilarity per sample)
sample_beta_diversity <- data.frame(
  Sample = rownames(community_matrix),
  Bray_Curtis_Mean = apply(bray_matrix, 1, function(x) mean(x[x > 0])),  # Exclude self-comparison (0)
  Jaccard_Mean = apply(jaccard_matrix, 1, function(x) mean(x[x > 0]))
)

# Prepare beta diversity for sample-wise plotting
beta_long <- sample_beta_diversity %>%
  pivot_longer(cols = c(Bray_Curtis_Mean, Jaccard_Mean), 
               names_to = "Index", 
               values_to = "Value") %>%
  mutate(Index = gsub("_Mean", "", Index),
         Index = gsub("_", "-", Index))

# =============================================================================
# Create Sample-wise Publication-Ready Plots
# =============================================================================

# Set theme for publication
pub_theme <- theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black")
  )

# Sample-wise Alpha Diversity Plot
alpha_sample_plot <- ggplot(alpha_long, aes(x = Sample, y = Value, fill = Index)) +
  geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
  geom_text(aes(label = round(Value, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 2.5, angle = 0) +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  facet_wrap(~Index, scales = "free_y", ncol = 2) +
  labs(
    title = "Alpha Diversity Indices by Sample",
    x = "Sample ID",
    y = "Diversity Index Value",
    #subtitle = paste(" ", length(unique(alpha_long$Sample)), " "),
    fill = "Diversity Index"
  ) +
  pub_theme +
  theme(legend.position = "bottom")

plot(alpha_sample_plot)
######

# Sample-wise Shannon Diversity (standalone plot)
shannon_sample_plot <- ggplot(alpha_diversity, aes(x = reorder(Sample, Shannon), y = Shannon, fill = Shannon)) +
  geom_col(alpha = 0.8, color = "black", size = 0.5) +
  geom_text(aes(label = round(Shannon, 3)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#3498db", high = "#e74c3c", name = "Shannon\nIndex") +
  labs(
    title = " ",# Shannon Diversity Index by Sample
    x = "Studies (Bovine Mastitis)",
    y = "Shannon Diversity Index (H')",
    subtitle = " "# whithin-sample species diversity
  ) +
  pub_theme +
  theme(legend.position = "right")

plot(shannon_sample_plot)

##############################################################################
# Sample-wise Beta Diversity Plot (Mean dissimilarity per sample)
beta_sample_plot <- ggplot(beta_long, aes(x = Sample, y = Value, fill = Index)) +
  geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
  geom_text(aes(label = round(Value, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 2.8, angle = 0) +
  scale_fill_manual(values = c("Bray-Curtis" = "#E31A1C", "Jaccard" = "#1F78B4")) +
  labs(
    title = "  ",
    x = "Studies",
    y = "Mean Dissimilarity",
    subtitle = "", #Average dissimilarity of each sample compared to all others
    fill = "Dissimilarity\nIndex"
  ) +
  pub_theme +
  theme(legend.position = "bottom")

plot(beta_sample_plot)

# Richness vs Shannon Plot
richness_shannon_plot <- ggplot(alpha_diversity, aes(x = Richness, y = Shannon)) +
  geom_point(aes(size = Simpson, color = Sample), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 10) +
  scale_size_continuous(range = c(3, 8), name = "Simpson\nIndex") +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "", #Species Richness vs Shannon Diversity
    x = "Species Richness (Number of Taxa)",
    y = "Shannon Diversity Index (H')",
    subtitle = "", #Relationship between richness and diversity (point size = Simpson index)
    color = "Studies"
  ) +
  pub_theme +
  theme(legend.position = "right")

plot(richness_shannon_plot)
#####
# Evenness by Sample Plot
evenness_plot <- ggplot(alpha_diversity, aes(x = reorder(Sample, Evenness), y = Evenness, fill = Evenness)) +
  geom_col(alpha = 0.8, color = "black", size = 0.5) +
  geom_text(aes(label = round(Evenness, 3)), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#e74c3c", mid = "#f39c12", high = "#2ecc71", 
                       midpoint = 0.5, name = "Evenness") +
  labs(
    title = "", #Species Evenness by Sample
    x = "Studies",#(ordered by evenness)
    y = "Pielou's Evenness Index (J')",
    subtitle = "" #How evenly species are distributed within each sample
  ) +
  pub_theme +
  theme(legend.position = "right")
plot(evenness_plot)
############
# Create output directory
if (!dir.exists("diversity_plots")) {
  dir.create("diversity_plots")
}

# Save sample-wise plots
ggsave("diversity_plots/Alpha_Diversity_by_Samplekk.png", alpha_sample_plot, 
       width = 14, height = 10, dpi = 300, bg = "white")

ggsave("diversity_plots/Shannon_Diversity_by_Sample.png", shannon_sample_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave("diversity_plots/Beta_Diversity_by_Sample.png", beta_sample_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")


ggsave("diversity_plots/Richness_vs_Shannon.png", richness_shannon_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave("diversity_plots/Evenness_by_Sample.png", evenness_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

##########################################################################
# Save PDF versions for publication
ggsave("diversity_plots/Alpha_Diversity_by_Sample.pdf", alpha_sample_plot, 
       width = 14, height = 10, device = "pdf")

ggsave("diversity_plots/Shannon_Diversity_by_Sample.pdf", shannon_sample_plot, 
       width = 12, height = 8, device = "pdf")

ggsave("diversity_plots/Beta_Diversity_by_Sample.pdf", beta_sample_plot, 
       width = 12, height = 8, device = "pdf")

ggsave("diversity_plots/Alpha_Diversity_Comparison.pdf", alpha_comparison_plot, 
       width = 14, height = 8, device = "pdf")

ggsave("diversity_plots/Richness_vs_Shannon.pdf", richness_shannon_plot, 
       width = 12, height = 8, device = "pdf")

ggsave("diversity_plots/Evenness_by_Sample.pdf", evenness_plot, 
       width = 12, height = 8, device = "pdf")

