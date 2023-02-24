
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
#BiocManager::install("GenomicRanges", force = T)
getwd()
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
rm(list = ls())
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/")

hub_data <- read.delim("Basic10kbProteinCodingContactHubsSimple120_promo.txt")

sCREs <- hub_data %>%
  filter(cCRE =="FALSE") 


#sCREs chr distribution

chr_count <- sCREs %>%
  group_by(Chr) %>%
  summarise(n_sCREs = n()) ###why only 20 chromosomes?


ggplot(data = sCREs, aes(x = Chr, group = Chr)) +
  geom_bar()


###exploration of gene function

gene_scre <-  sCREs %>%
  dplyr::select(CisBin,Genes) 


#gene_list <- (separate(data = gene_ccre, col = Genes, into = c(as.character(seq(1,41,1))),","))
write.csv(gene_scre, file = "list_gene_scres.csv")
save(file = "gene_list_scres.rda", gene_scre)


gene_list <- gene_list %>%
  select(-CisBin)
gene_list_long <- gather(gene_list) 
gene_list_long <- na.omit(gene_list_long)
gene_number <- c(seq(6000,24000, 6000))
# 
# gene_list_one = gene_list_long[1:gene_number[1],]
# gene_list_two = gene_list_long[gene_number[1]:gene_number[2],]
# gene_list_three = gene_list_long[gene_number[2]:gene_number[3],]
# gene_list_four = gene_list_long[gene_number[3]:gene_number[4],]
# write.csv(gene_list_long_1, file = "list_genes_one.csv")
# write.csv(gene_list_two, file = "list_genes_two.csv")
# write.csv(gene_list_three, file = "list_genes_three.csv")
# write.csv(gene_list_four, file = "list_genes_four.csv")

functions <- read.delim("tr_646566BCCBFA1674162501012.txt")
#my.dnastring <- as.character(Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, "chr1", 3000000, 3000100))
mouse.genome <- BSgenome.Mmusculus.UCSC.mm39


#test<- getSeq(BSgenome.Mmusculus.UCSC.mm39, 'chr1', start=16650000, end=16660000)
sequence_list <- list()
for (i in 1:nrow(sCREs)) {
   sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm39, paste(sCREs[i,6]), start=sCREs[i,7], end=sCREs[i,8])
   sequence_list[[i]] <-sequence
}


save(file = "sequence_list.rda", sequence_list)


load("sequence_list.rda")
sequence_df <- data.frame()
for (i in 1:length(sequence_list)) { 
  #extracts each sequence
  one_sequence <- sequence_list[[i]]
  #convers each sequence to a df and adds it to a df
  sequence_df[i,1] <-data.frame(toString(one_sequence))
}

sCREs <- cbind(sCREs,sequence_df)
sCREs <- sCREs %>%
  rename("toString.one_sequence." = "sequence")

sCREs_export <- sCREs %>%
  select(CisBin,sequence)
D <- do.call(rbind, lapply(seq(nrow(sCREs_export)), function(i) t(sCREs_export[i, ])))


bed_format <- sCREs %>%
  select(CisBin, Chr, Start, End) %>%
  rename("Start" = "chromStart",
         "End" = "chromEnd",
         "Chr" = "chrom",
         "CisBin" = "name")

write.csv(bed_format, file = "file.csv")
GenomicRanges::makeGRangesFromDataFrame




library(memes)
runStreme(D, control ="shuffle")
