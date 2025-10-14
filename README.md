# Brassinosteroid_XTH_bean
Data and R scripts for: Brassinosteroid treatment reveals the importance of xyloglucan transglucosylase/hydrolase (XTH) genes in growth habit determination of twining common bean vines.

Note that the three treatments are abreviated as BL for brassinolide-treated (a brassinosteroid), BZ for brassinazole treated (a brassinosteroid inhibitor), and ET for the ethanol-only "mock" treatment (lanolin with the equivalent amount of ethanol used to dissove brassinolide and brassinazole). 

The R codes used to generate visuals are listed below with the Figure from the paper.

When running codes from the RNAseq data, it is important to run them in sequence:

First, convert the .BAM files into gene_counts and metadata. .BAM files can be downloaded from (website) and the resulting gene_counts.xlsx and metadata.xlsx are available for download to run the next steps if you choose to analyze the data from there.
Next, 
Third, run the DEG and volcano plots code which will provide the dds necessary for the next step (heatmap)
Fourth,  XTH Gene Expression Heatmap.


Guide to Scripts by Figure

Figure 1: 
-Internode 4 Elongation

Figure 2: 
-G-Fiber Thickness
-G-fiber Length

Figure 3: Venn Diagram and GO terms

Figure 4: 
-Volcano Plots
-Xyloglucan

Figure 5:
-XTH Heatmap

Figure 6: (no codes used)
