#### STEP 1: Install and Load Packages
packages <- c("VennDiagram", "ggplot2", "DESeq2", "Rsamtools", "GenomicAlignments", "tximport", "edgeR", "Rsubread")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")  # Ensure you're on the correct Bioconductor version

BiocManager::install(c("Rsamtools", "GenomicAlignments", "tximport", "edgeR", "Rsubread"), ask=FALSE)
BiocManager::install("edgeR")



#### STEP 2: Create Sample Table
bam_dir <- "C:/Users/lenam/OneDrive/Documents/Beans/RNAseq Results/02.Bam_Oct2024"
bam_files <- list.files(bam_dir, pattern = "*.bam$", full.names = TRUE)

# --- FIXED extraction of time_point and treatment ---
# Example expected filename pattern: "A_BL_1.bam" or "B_ET_2.bam"
sample_table <- data.frame(
  file = bam_files,
  time_point = sub("^([A-Z])_.*$", "\\1", basename(bam_files)),             # Extract the first letter before first underscore
  treatment = sub("^[A-Z]_([A-Z]+)_\\d+\\.bam$", "\\1", basename(bam_files)), # Extract treatment (middle part)
  stringsAsFactors = FALSE
)

#### STEP 3: Count Reads Using FeatureCounts
gtf_file <- "C:/Users/lenam/OneDrive/Documents/Beans/RNAseq Results/04.Ref/genome.gtf/genome.gtf"

counts <- featureCounts(files = bam_files, 
                        annot.ext = gtf_file,        
                        isGTFAnnotationFile = TRUE,  
                        GTF.featureType = "exon",    
                        GTF.attrType = "gene_id",    
                        useMetaFeatures = TRUE,      
                        isPairedEnd = TRUE,          
                        nthreads = 4)

gene_counts <- counts$counts  # Directly assign

# Clean row names
rownames(gene_counts) <- sub("gene_PHAVU_", "", rownames(gene_counts))

# Clean column names
colnames(gene_counts) <- sub("\\.bam$", "", basename(colnames(gene_counts)))


# Prepare metadata with corrected row names
metadata <- data.frame(
  row.names = colnames(gene_counts),
  time_point = sample_table$time_point,
  treatment = sample_table$treatment
)

## At this point: gene_counts: Rows = Gene IDs (with "gene_PHAVU_" removed) 
## Columns = Sample names (cleaned), and metadata aligned properly.



# Filtering gene_counts! Removing genes with < 10 counts per million (CPM)

library(edgeR)

# Convert gene_counts to a DGEList object
dge <- DGEList(counts = gene_counts)

# Calculate CPM values
cpm_values <- cpm(dge)

# Group Samples by Treatment and Timepoint
sample_groups <- paste(metadata$treatment, metadata$time_point, sep = "_")

# Get unique groups
unique_groups <- unique(sample_groups)
print(unique_groups)

# Apply Filtering Per Group
keep_genes_per_group <- sapply(unique_groups, function(group) {
  sample_indices <- which(sample_groups == group)
  mat <- cpm_values[, sample_indices, drop = FALSE]
  matrixStats::rowMins(mat) >= 10
})

keep_genes <- rowSums(keep_genes_per_group) > 0

# Apply filtering
filtered_gene_counts <- gene_counts[keep_genes, ]

# Check dimensions before and after filtering
dim(gene_counts)
dim(filtered_gene_counts)

# Rewrite gene_counts with filtered data
gene_counts <- filtered_gene_counts

# Confirm matrix type
class(gene_counts)



# Install and load writexl if needed
if(!require(writexl)) install.packages("writexl")
library(writexl)

# -----------------------------
# Export filtered gene counts to Excel
# -----------------------------
gene_counts_df <- data.frame(
  GeneID = rownames(gene_counts),
  gene_counts,
  check.names = FALSE
)
write_xlsx(gene_counts_df, "gene_counts.xlsx")

# -----------------------------
# Export metadata to Excel
# -----------------------------
metadata_df <- data.frame(
  SampleID = rownames(metadata),
  metadata,
  check.names = FALSE
)
write_xlsx(metadata_df, "metadata.xlsx")
