# the script reads in 4 CSV files and combines the data into one dataframe. 
# Then it performs a gene ontology (GO) enrichment analysis on the entrez IDs in that dataframe. 
# Finally, it creates a ggplot of the enrichment results showing the 
# Rich Factor (based on a ratio of the number of entrez IDs in a particular GO term to the total number of entrez IDs) 
# the x-axis and the reordered Description column from the enrichment results on the y-axis. 
# The plot uses color and size to indicate the p-value and count, respectively.
#Code was commented using ChatGPT
#Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation 2, 100141 (2021).

library(ggplot2) # load ggplot2 library for data visualization
library(reshape2) # load reshape2 library for reshaping data
library(tidyverse) # load tidyverse library for data manipulation and visualization
library(gridExtra) # load gridExtra library for creating complex grid-based layouts
library(data.table) # load data.table library for manipulating large data sets
library(ggpubr) # load ggpubr library for creating publication-ready plots
library(ggridges) # load ggridges library for creating density ridgeline plots
library(BSgenome.Mmusculus.UCSC.mm39) # load BSgenome.Mmusculus.UCSC.mm39 library for mouse genome data
library(seqinr) # load seqinr library for bioinformatics
library(Biostrings) # load Biostrings library for handling DNA sequences
library(GenomicRanges) # load GenomicRanges library for handling genomic intervals
library(GOSemSim) # load GOSemSim library for calculating semantic similarity between Gene Ontology terms
library(clusterProfiler) # load clusterProfiler library for gene set enrichment analysis
#BiocManager::install("disgenet2r")
getwd() # get the current working directory
rm(list = ls()) # clear the workspace of all objects in the current environment
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/") # set the working directory

gene_list_1 <- read.csv("list_genes_one.csv") # read file list_genes_one.csv and store it in gene_list_1
gene_list_2 <- read.csv("list_genes_two.csv") # read file list_genes_two.csv and store it in gene_list_2
gene_list_3 <- read.csv("list_genes_three.csv") # read file list_genes_three.csv and store it in gene_list_3
gene_list_4 <- read.csv("list_genes_four.csv") # read file list_genes_four.csv and store it in gene_list_4
gene_list <- (data.frame(rbind(gene_list_1,gene_list_2,gene_list_3,gene_list_4))) # combine data from gene_list_1, gene_list_2, gene_list_3, gene_list_4 and store it in gene_list
gene_list <- gene_list %>%
  dplyr::select(value) # select the value column of the gene_list

library(org.Mm.eg.db) # load org.Mm.eg.db library
mm <- org.Mm.eg.db # assign mm to org.Mm.eg.db
my.symbols <- gene_list$value # assign my.symbols to the value column of gene_list_1
gene_entrez <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") # select the ENTREZID and SYMBOL columns from mm for the symbols in my.symbols and store it in gene_entrez
ego <- enrichGO(gene_entrez$ENTREZID, OrgDb= "org.Mm.eg.db", ont="BP", readable=TRUE) # perform gene ontology (GO) enrichment analysis on the ENTREZID values in gene_entrez using the org.Mm.eg.db database and store the result in ego
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min) # simplify the ego object by selecting the lowest p-value for each GO term
ego3 <- mutate(ego, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio))) # create a new column "richFactor" in ego by dividing the "Count" column by the background ratio

ggplot(ego3, showCategory = 10, aes(richFactor, fct_reorder(Description, richFactor))) + # create a ggplot of ego3 with the x-axis as "richFactor" and the y-axis as the reordered "Description" column
  geom_segment(aes(xend=0, yend = Description)) + # add a segment on the plot connecting the x-axis to each point
  geom_point(aes(color=p.adjust, size = Count)) + # add points to the plot colored by p-value and sized by count
  scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) + # set the color scale for the points
  scale_size_continuous(range=c(2, 10)) + # set the size scale for the points
  xlab("Rich Factor") + # label the x-axis
  ggtitle("Biological Processes/cCREs") # add a title to the plot

#####Gene/disease relationship
#BiocManager::install( "biomaRt" )

save(file = "ego2_function.rda", ego2)
save(file = "ego3_function.rda", ego3)
