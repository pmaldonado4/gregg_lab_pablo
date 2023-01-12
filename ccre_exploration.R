
library(ggplot2)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(ggpubr)
library(ggridges)

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

gene_list_long <- gather(gene_list,CisBin)


