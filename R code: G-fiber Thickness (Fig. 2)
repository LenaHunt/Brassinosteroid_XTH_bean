# Load the necessary libraries
library(ggplot2)
library(dplyr)

# Define colors for each Treatment
colors <- c("ET" = "darkturquoise", "BL" = "orchid", "BZ" = "burlywood4")
dot_colors <- c("ET" = "turquoise", "BL" = "orchid2", "BZ" = "burlywood3")

# Filter out outliers (values above 9)
G_fiber_Thickness <- G_fiber_Thickness %>%
  filter(Length_um <= 9)

# Reorder the Treatment factor (keep ET, BL, BZ as-is)
G_fiber_Thickness <- G_fiber_Thickness %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("ET", "BL", "BZ")))

# Plotting
ggplot(G_fiber_Thickness, aes(x = Treatment, y = Length_um, fill = Treatment)) +
  # Add the violin plot
  geom_violin(outlier.shape = NA, width = 0.5, color = "black") +
  # Add jittered dots
  geom_jitter(aes(color = Treatment), width = 0.2, size = 2, alpha = 0.7) +
  # Add boxplot overlay
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  # Customizations
  labs(title = "G-fiber Thickness", y = "Length in µm") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "top"                  # Legend at the top
  ) +
  scale_fill_manual(values = colors) +       # Custom fill colors
  scale_color_manual(values = dot_colors)    # Custom colors for dots


### STATS

# Perform one-way ANOVA
anova_result <- aov(Length_um ~ Treatment, data = G_fiber_Thickness)
anova_summary <- summary(anova_result)
anova_p_value <- anova_summary[[1]]$`Pr(>F)`[1]

# Perform *unpaired* t-tests between each pair of groups
t_test_ET_BL <- t.test(G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "ET"],
                       G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "BL"],
                       paired = FALSE)

t_test_ET_BZ <- t.test(G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "ET"],
                       G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "BZ"],
                       paired = FALSE)

t_test_BL_BZ <- t.test(G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "BL"],
                       G_fiber_Thickness$Length_um[G_fiber_Thickness$Treatment == "BZ"],
                       paired = FALSE)

# Extract p-values
p_value_ET_BL <- t_test_ET_BL$p.value
p_value_ET_BZ <- t_test_ET_BZ$p.value
p_value_BL_BZ <- t_test_BL_BZ$p.value

# Print the p-values
cat("ANOVA p-value:", anova_p_value, "\n")
cat("Unpaired t-test ET vs BL p-value:", p_value_ET_BL, "\n")
cat("Unpaired t-test ET vs BZ p-value:", p_value_ET_BZ, "\n")
cat("Unpaired t-test BL vs BZ p-value:", p_value_BL_BZ, "\n")

