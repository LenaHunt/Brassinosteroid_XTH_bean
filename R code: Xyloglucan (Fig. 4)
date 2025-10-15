# Load libraries
library(ggplot2)

# Custom colors for the treatments (update "ET" instead of "mock")
custom_colors <- c("ET" = "darkturquoise", "BL" = "orchid", "BZ" = "burlywood4")
dot_colors <- c("ET" = "turquoise", "BL" = "orchid2", "BZ" = "burlywood3")

# Reorder Treatment factor without renaming ET
Xyloglucan$Treatment <- factor(Xyloglucan$Treatment, levels = c("ET", "BL", "BZ"))

# Violin plot with jittered dots
ggplot(Xyloglucan, aes(x = Treatment, y = Xyloglucan, fill = Treatment)) +
  geom_violin(outlier.shape = NA, width = 0.5, color = "black") +  # Violin plot
  geom_jitter(aes(color = Treatment), width = 0.1, size = 2, alpha = 0.6) +  # Jittered dots
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +  # Box plot overlay
  labs(title = "Xyloglucan by Treatment", y = "Xyloglucan speckle area") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = dot_colors)

# Print the number of entries for each Treatment
table(Xyloglucan$Treatment)

# Perform ANOVA
anova_model <- aov(Xyloglucan ~ Treatment, data = Xyloglucan)
summary(anova_model)

# Perform t-tests between each pair of treatments with Bonferroni correction
pairwise_t_test <- pairwise.t.test(Xyloglucan$Xyloglucan, Xyloglucan$Treatment, p.adjust.method = "bonferroni")
pairwise_t_test
