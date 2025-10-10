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

sample_table <- data.frame(
  file = bam_files,
  time_point = substr(basename(bam_files), 1, 1),
  treatment = ifelse(grepl("_C\\d+\\.bam$", basename(bam_files)), "C", 
                     sub(".*_([A-Z]+)_\\d+\\.bam$", "\\1", basename(bam_files))),
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

## At this point: gene_counts: Rows = Gene ID number (with "gene_PHAVU_" removed) 
## and Columns = Sample names (cleaned to only contain sample identifiers)



# Filtering gene_counts! Removing genes with < 10 counts per million (CPM)

library(edgeR)

# Convert gene_counts to a DGEList object
dge <- DGEList(counts = gene_counts)

# Calculate CPM values
cpm_values <- cpm(dge)

# Group Samples by Treatment and Timepoint

# Extract groups based on metadata
sample_groups <- paste(metadata$treatment, metadata$time_point, sep = "_")

# Get unique groups
unique_groups <- unique(sample_groups)

print(unique_groups)

# Apply Filtering Per Group

# Create a logical matrix to store whether each gene meets the 10 CPM threshold per group
keep_genes_per_group <- sapply(unique_groups, function(group) {
  # Get sample indices for this group
  sample_indices <- which(sample_groups == group)
  
  # Subset the columns as a matrix, even if there's only one sample
  mat <- cpm_values[, sample_indices, drop = FALSE]
  
  # Check if the gene has at least 10 CPM in **all** samples of this group
  matrixStats::rowMins(mat) >= 10
})

# Keep genes that pass the filter in **at least one group**
keep_genes <- rowSums(keep_genes_per_group) > 0



# Apply filtering
filtered_gene_counts <- gene_counts[keep_genes, ]

# Check dimensions of the filtered dataset (12586 genes in 67 samples)
dim(filtered_gene_counts)

dim(gene_counts)  # Before filtering
dim(filtered_gene_counts)  # After filtering

#rewrite the old gene_counts
gene_counts <- filtered_gene_counts

#Download to computer
class(gene_counts)



# Install and load writexl if needed
if(!require(writexl)) install.packages("writexl")
library(writexl)

# -----------------------------
# Export filtered gene counts to Excel
# -----------------------------
# Convert matrix to data frame and add GeneID as first column
gene_counts_df <- data.frame(
  GeneID = rownames(gene_counts),
  gene_counts,
  check.names = FALSE
)

# Write to Excel
write_xlsx(gene_counts_df, "gene_counts.xlsx")




