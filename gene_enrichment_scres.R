
# It appears that the code performs a Gene Ontology enrichment analysis on a list of genes, 
# then simplifies the results and saves them to a file. 
# Finally, it creates a plot of the results using the ggplot2 library. 
# The plot shows the "rich factor" of the genes in relation to their corresponding biological process, 
# and the size and color of the points on the plot correspond to the 
# "count" and "p.adjust" values of the genes, respectively.
#Commented using ChatGPT
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(ggpubr)
library(ggridges)
library(BSgenome.Mmusculus.UCSC.mm39)
library(seqinr)
library(Biostrings)
library(GenomicRanges)
library(GOSemSim)
library(clusterProfiler)
library(GOSemSim)

# Install GenomicRanges if not already installed
#BiocManager::install("GenomicRanges", force = T)

# Set the working directory
getwd()
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
rm(list = ls())
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/")

# Load the gene list 
load("gene_list_scres.rda")
gene_list <- gene_scre

#Select the Genes column of the gene list
gene_list <- gene_list %>%
  dplyr::select(Genes)

# Load the org.Mm.eg.db library and select entrez ID and symbol for the genes in the list
library(org.Mm.eg.db)
mm <- org.Mm.eg.db
my.symbols <- gene_list$Genes
gene_entrez <- select(mm, 
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")

# Perform Gene Ontology enrichment analysis
ego <- enrichGO(gene_entrez$ENTREZID, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)

# Simplify the results of the analysis
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

# Create a new variable 'richFactor'
ego3 <- mutate(ego, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))

# Save the results to a file
save(file="ego3.rda", ego3)

# Create a plot of the results
ggplot(ego3, showCategory = 10, aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) + 
  scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) + 
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor") + 
  ggtitle("Biological Processes/sCREs")
