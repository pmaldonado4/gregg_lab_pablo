#Step 1 get DNA sequences  using GRanges
library(GenomicRanges)
library(org.Mm.eg.db) # load org.Mm.eg.db library
library(memes)
library(tidyverse)
library(universalmotif)
check_meme_install()
mouse.genome <- BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39
rm(list = ls()) # clear the workspace of all objects in the current environment
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Active Projects/Gregg Lab/cCREs Project/") # set the working directory
# You can manually input a path to meme_path
# If no meme/bin is detected, will return a red X
check_meme_install(meme_path = "/opt/local/bin")
options(meme_bin = "/opt/local/bin") 
hub_data <- read.delim("Basic10kbProteinCodingContactHubsSimple120_promo.txt") %>%
  mutate(range = paste0(Chr,sep = ":",Start, sep = "-", End))
test<- get_sequence(hub_data[1:100,10], mouse.genome)
dreme_results <- runDreme(test, control = "shuffle")

dreme_results %>% 
  to_list() %>% 
  view_motifs()
