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
library(org.Mm.eg.db) # load org.Mm.eg.db library
#BiocManager::install("disgenet2r")
getwd() # get the current working directory
rm(list = ls()) # clear the workspace of all objects in the current environment
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/") # set the working directory
mm <- org.Mm.eg.db # assign mm to org.Mm.eg.db

gene_list_1 <- read.csv("list_genes_one.csv") # read file list_genes_one.csv and store it in gene_list_1
gene_list_2 <- read.csv("list_genes_two.csv") # read file list_genes_two.csv and store it in gene_list_2
gene_list_3 <- read.csv("list_genes_three.csv") # read file list_genes_three.csv and store it in gene_list_3
gene_list_4 <- read.csv("list_genes_four.csv") # read file list_genes_four.csv and store it in gene_list_4
gene_list <- (data.frame(rbind(gene_list_1,gene_list_2,gene_list_3,gene_list_4))) # combine data from gene_list_1, gene_list_2, gene_list_3, gene_list_4 and store it in gene_list
gene_list <- gene_list %>% 
  dplyr::select(value) # select the value column of the gene_list

#load all scre and ccre gene data
hub_data <- read.delim("Basic10kbProteinCodingContactHubsSimple120_promo.txt")
gene_all <-  hub_data %>%
  dplyr::select(CisBin,Genes) 
genes_brain <- read.csv("brain_expressed_genes.csv",row.names=NULL)
gene_list_all <- genes_brain$CD.ID..Child.



gene_list_long_all<- unique(gene_list_all)

gene_entrez_all <- AnnotationDbi::select(mm, keys = genes_brain$CD.ID..Child., columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") 


my.symbols <- gene_list$value # assign my.symbols to the value column of gene_list_1

gene_entrez <- AnnotationDbi::select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") # select the ENTREZID and SYMBOL columns from mm for the symbols in my.symbols and store it in gene_entrez
gene_entrez <- gene_entrez %>%
  distinct()
ego <- enrichGO(gene_entrez$ENTREZID, OrgDb= "org.Mm.eg.db", ont="BP", readable=TRUE, pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05) # perform gene ontology (GO) enrichment analysis on the ENTREZID values in gene_entrez using the org.Mm.eg.db database and store the result in ego
#drop_it <- filterGO(ego)

ego2 <- simplify(ego, cutoff=0.3, by="p.adjust", select_fun=min) # simplify the ego object by selecting the lowest p-value for each GO term

ego3 <- mutate(ego2, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio))) # create a new column "richFactor" in ego by dividing the "Count" column by the background ratio


ggplot(ego3, showCategory = 20, aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2",
                    "#7e62a3"),
                    trans = "log10",
  guide=guide_colorbar(reverse=TRUE,
                       order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("Biological Processes") +
  theme_classic()
  
  
  
save(file = "ego2_function.rda", ego2)

save(file = "ego3_function.rda", ego3)


ggsave("enrichment_all_genes_top_20_ccres_brain_universe.pdf", width = 10,height = 10, dpi = 700)


#####Gene/disease relationship
#BiocManager::install( "biomaRt" )

save(file = "ego2_function.rda", ego2)
save(file = "ego3_function.rda", ego3)
