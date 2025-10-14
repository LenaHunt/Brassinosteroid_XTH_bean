# Brassinosteroid_XTH_bean
Data and R scripts for: Brassinosteroid treatment reveals the importance of xyloglucan transglucosylase/hydrolase (XTH) genes in growth habit determination of twining common bean vines.

Note that the three treatments are abreviated as BL for brassinolide-treated (a brassinosteroid), BZ for brassinazole treated (a brassinosteroid inhibitor), and ET for the ethanol-only "mock" treatment (lanolin with the equivalent amount of ethanol used to dissove brassinolide and brassinazole). 

The R codes used to generate visuals are listed below with the Figure from the paper.

When running codes from the RNAseq data, it is important to first convert the .BAM files into gene_counts and metadata which can be used in the next steps. The order in which the codes should be run are: 1.) BAM_fules_to_gene_counts, 2.) DEG and volcano plots, 3.) XTH Gene Expression Heatmap


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
