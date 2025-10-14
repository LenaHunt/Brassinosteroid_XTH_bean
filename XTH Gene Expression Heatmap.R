# Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# 1. VST transform
vst_data <- vst(dds, blind = FALSE)
vst_counts <- assay(vst_data)


# 2. Define your genes of interest
XTH_genes_of_interest <- c(
  "004G168700", "002G311200", "009G233200", "003G137600", "003G147500",
  "003G147600", "003G147400", "003G143332", "003G147700", "003G147300",
  "002G244000", "001G265300", "006G033800", "008G129100", "008G129900",
  "007G148500", "001G179900", "005G130900", "011G085200", "003G052400",
  "002G136200", "002G182000", "003G014400", "005G111300", "011G107000",
  "003G231900", "002G008400", "001G098200", "001G098300", "004G062400",
  "004G111800", "007G052800", "007G048300"
)


# Filter for selected genes
selected_counts_XTH <- vst_counts[rownames(vst_counts) %in% XTH_genes_of_interest, , drop = FALSE]

# Add missing genes as zero rows
missing_genes <- setdiff(XTH_genes_of_interest, rownames(selected_counts_XTH))
if(length(missing_genes) > 0){
  zero_matrix <- matrix(0, nrow = length(missing_genes), ncol = ncol(vst_counts))
  rownames(zero_matrix) <- missing_genes
  colnames(zero_matrix) <- colnames(vst_counts)
  selected_counts_XTH <- rbind(selected_counts_XTH, zero_matrix)
}

# Define the full custom order
custom_order <- XTH_genes_of_interest  # already in desired order

# Reorder rows according to custom order
selected_counts_XTH <- selected_counts_XTH[custom_order, , drop = FALSE]


# 3. Aggregate samples by treatment group (trimmed mean 10%)
sample_groups <- gsub("_[0-9]+$", "", colnames(selected_counts_XTH))  # e.g., "A_ET"
aggregated_counts <- sapply(unique(sample_groups), function(group) {
  apply(selected_counts_XTH[, sample_groups == group, drop = FALSE], 1, mean, trim = 0.1)
})
aggregated_counts <- as.matrix(aggregated_counts)

# Replace any invalid values with 0
aggregated_counts[!is.finite(aggregated_counts)] <- 0


# 4. Define heatmap colors
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Desired column order
desired_order <- c("A_ET","A_BL","A_BZ","B_ET","B_BL","B_BZ")
aggregated_counts <- aggregated_counts[, desired_order, drop = FALSE]


# 5. Generate heatmap
pheatmap(
  aggregated_counts,
  color = heatmap_colors,
  main = "XTH Gene Expression Heatmap (VST, Custom Order with Zero Rows)",
  cluster_rows = FALSE,   # preserve custom order
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row"
)
