# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Load your dataset
data <- G_fiber_Length

# Reshape your data
long_data <- pivot_longer(data, cols = c("BL", "BZ", "ET"), names_to = "Group", values_to = "Values")

# Calculate mean and standard error
summary_data <- long_data %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Values),
    SE = sd(Values) / sqrt(n())
  )


# Reorder the levels of the Group variable
long_data$Group <- factor(long_data$Group, levels = c("ET", "BL", "BZ"))

# Plotting
ggplot(long_data, aes(x = Group, y = Values, fill = Group)) +
  # Add the violin plot
  geom_violin(outlier.shape = NA, width = 0.5, color = "black") +
  # Add jittered dots
  geom_jitter(aes(color = Group), width = 0.2, size = 2, alpha = 0.7) +
  # Add boxplot overlay
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  # Customizations
  labs(title = "G-fiber length", y = "Length in mm") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "top"  # Position the legend at the top
  ) +
  scale_fill_manual(values = c("ET" = "darkturquoise", "BL" = "orchid", "BZ" = "burlywood4")) +  # Custom colors for boxplot fill
  scale_color_manual(values = c("ET" = "turquoise", "BL" = "orchid2", "BZ" = "burlywood3"))  # Custom colors for dots


### MOCK ###


# Plotting with updated label for ET
ggplot(long_data, aes(x = Group, y = Values, fill = Group)) +
  # Add the violin plot
  geom_violin(outlier.shape = NA, width = 0.5, color = "black") +
  # Add jittered dots
  geom_jitter(aes(color = Group), width = 0.2, size = 2, alpha = 0.7) +
  # Add boxplot overlay
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  # Customizations
  labs(title = "G-fiber length", y = "Length in mm") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "top"  # Position the legend at the top
  ) +
  scale_fill_manual(values = c("ET" = "darkturquoise", "BL" = "orchid", "BZ" = "burlywood4")) +  # Custom colors for boxplot fill
  scale_color_manual(values = c("ET" = "turquoise", "BL" = "orchid2", "BZ" = "burlywood3")) +  # Custom colors for dots
  scale_x_discrete(labels = c("ET" = "mock", "BL" = "BL", "BZ" = "BZ"))  # Change label for ET



####STATS

# Perform ANOVA
anova_result <- aov(Values ~ Group, data = long_data)
anova_summary <- summary(anova_result)
anova_p_value <- anova_summary[[1]]$`Pr(>F)`[1]

# Perform paired t-tests between each pair of groups
t_test_BL_BZ <- t.test(long_data$Values[long_data$Group == "BL"], 
                       long_data$Values[long_data$Group == "BZ"], 
                       paired = TRUE)
t_test_BL_ET <- t.test(long_data$Values[long_data$Group == "BL"], 
                       long_data$Values[long_data$Group == "ET"], 
                       paired = TRUE)
t_test_BZ_ET <- t.test(long_data$Values[long_data$Group == "BZ"], 
                       long_data$Values[long_data$Group == "ET"], 
                       paired = TRUE)

# Extract p-values
p_value_BL_BZ <- t_test_BL_BZ$p.value
p_value_BL_ET <- t_test_BL_ET$p.value
p_value_BZ_ET <- t_test_BZ_ET$p.value

# Print the p-values
cat("ANOVA p-value:", anova_p_value, "\n")
cat("Paired t-test BL vs BZ p-value:", p_value_BL_BZ, "\n")
cat("Paired t-test BL vs ET p-value:", p_value_BL_ET, "\n")
cat("Paired t-test BZ vs ET p-value:", p_value_BZ_ET, "\n")
