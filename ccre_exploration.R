
library(ggplot2)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(ggpubr)
library(ggridges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
getwd()

rm(list = ls())
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/")

hub_data <- read.delim("Basic10kbProteinCodingContactHubsSimple120_promo.txt")

cCREs <- hub_data %>%
  filter(cCRE =="TRUE") 


#ccres chr distribution

chr_count <- cCREs %>%
  group_by(Chr) %>%
  summarise(n_ccres = n()) ###why only 20 chromosomes?


ggplot(data = cCREs, aes(x = Chr, group = Chr)) +
  geom_bar()


###exploration of gene function

gene_ccre <-  cCREs %>%
  select(CisBin,Genes) 


gene_list <- (separate(data = gene_ccre, col = Genes, into = c(as.character(seq(1,41,1))),","))

gene_list_long <- gather(gene_list)
my.dnastring <- as.character(Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, "chr1", 3000000, 3000100))

#test<- getSeq(BSgenome.Mmusculus.UCSC.mm10, 'chr1', start=16650000, end=16660000)
# sequence_list <- list()
# for (i in 1:nrow(cCREs)) {
#   sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, paste(cCREs[i,6]), start=cCREs[i,7], end=cCREs[i,8])
#   sequence_list[[i]] <-sequence
# }


save(file = "sequence_list.rda", sequence_list)


load("sequence_list.rda")
sequence_df <- data.frame()
for (i in 1:length(sequence_list)) { 
  #extracts each sequence
  one_sequence <- sequence_list[[i]]
  #convers each sequence to a df and adds it to a df
  sequence_df[i,1] <-data.frame(toString(one_sequence))
}

cCREs <- cbind(cCREs,sequence_df)
cCREs <- cCREs %>%
  rename("toString.one_sequence." = "sequence")

cCREs_export <- cCREs %>%
  select(CisBin,sequence)
D <- do.call(rbind, lapply(seq(nrow(cCREs_export)), function(i) t(cCREs_export[i, ])))


bed_format <- cCREs %>%
  select(CisBin, Chr, Start, End) %>%
  rename("Start" = "chromStart",
         "End" = "chromEnd",
         "Chr" = "chrom",
         "CisBin" = "name")

write.csv(bed_format, file = "file.csv")
GenomicRanges::makeGRangesFromDataFrame