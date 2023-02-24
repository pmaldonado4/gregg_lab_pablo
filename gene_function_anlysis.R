BiocManager::install("GOSemSim")
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
#BiocManager::install("GenomicRanges", force = T)
getwd()
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
rm(list = ls())
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/")

list_one <- read.delim("tr_646566BCCBFA1674162501012.txt")
gene_list_1 <- read.csv("list_genes_one.csv")
gene_list_1 <- gene_list_1 %>%
  select(value)
keytypes(org.Dm.eg.db)
enrichGO(gene = gene_list_1,
         universe = names(gene_list),
         OrgDb = organism, 
         keyType = SYMBOL,
         readable = T,
         ont = "BP",
         pvalueCutoff = 0.05, 
         qvalueCutoff = 0.10)





list_two <- read.delim("tr_646566BCCBFA1674163541784.txt")
list_three <- read.delim("tr_646566BCCBFA1674163599607.txt")
list_four <- read.delim("tr_646566BCCBFA1674163649432.txt")
complete_list <- rbind(list_one,list_two,list_three,list_four)


complete_list <- complete_list %>%
  select(ID, Gene.Name,GOTERM_BP_DIRECT,UP_KW_BIOLOGICAL_PROCESS,UP_KW_DISEASE)



gene_function <- complete_list %>% 
  separate_rows(GOTERM_BP_DIRECT, sep=",") %>% 
  separate(GOTERM_BP_DIRECT, into=c("Function", "GOTERM_BP_DIRECT"),sep = "_", convert = TRUE) %>%
  select(-GOTERM_BP_DIRECT)
gene_function <- gene_function[-which(gene_function$Function == ""), ]
ggplot(data = gene_function) + 
geom_bar(aes(x = Function), fill = "blue") +
  xlab("Function") + ylab("Count")


length(unique(gene_function$Function))
gene_list_1 <- read.csv("list_genes_one.csv")

sim<-mgeneSim(gene_list_1,ont="MF",
              organism="mouse",measure="Wang")
