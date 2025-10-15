# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define custom colors
custom_colors <- c("ET" = "darkturquoise", "BL" = "orchid", "BZ" = "burlywood4")
custom_colors_points <- c("ET" = "turquoise1", "BL" = "orchid1", "BZ" = "burlywood2")

# Filter and prepare data
P5_combined <- Internode_4_Elongation %>%
  filter(Timepoint %in% c(0, 1, 2), Treatment %in% c("ET", "BL", "BZ"), Internode == "int4") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("ET", "BL", "BZ")),
    # Factor ordered by timepoint first, then treatment
    Treatment_Timepoint = factor(paste(Treatment, Timepoint, sep = "-"),
                                 levels = c(
                                   "ET-0", "BL-0", "BZ-0",
                                   "ET-1", "BL-1", "BZ-1",
                                   "ET-2", "BL-2", "BZ-2"
                                 ))
  )

# Create the plot
ggplot(P5_combined, aes(x = Treatment_Timepoint, y = Length, fill = Treatment)) +
  geom_violin(outlier.shape = NA, width = 0.8, color = "black") +
  geom_jitter(aes(color = Treatment), width = 0.2, size = 2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors_points) +
  labs(
    title = "Internode 4 - Length Across Timepoints",
    x = "Treatment-Timepoint",
    y = "Length (µm)",
    fill = "Treatment",
    color = "Treatment"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

##### STATS Kruskal-Wallis + pairwise Wilcoxon


# Function for non-parametric analysis per timepoint
analyze_nonparametric <- function(data, timepoint) {
  cat("\n\n#### Analysis for Timepoint", timepoint, "####\n")
  
  tp_data <- data %>% filter(Timepoint == timepoint)
  
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(Length ~ Treatment, data = tp_data)
  print(kruskal_result)
  
  if (kruskal_result$p.value < 0.05) {
    # Pairwise Wilcoxon test with Bonferroni correction
    wilcox_result <- pairwise.wilcox.test(tp_data$Length, tp_data$Treatment,
                                          p.adjust.method = "bonferroni",
                                          exact = FALSE)
    print(wilcox_result)
  } else {
    cat("No significant differences between groups.\n")
  }
}

# Run non-parametric analysis for each timepoint
for (tp in c(0, 1, 2)) {
  analyze_nonparametric(P5_combined, tp)
}
