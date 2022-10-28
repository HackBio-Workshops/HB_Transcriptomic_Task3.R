# DE Analysis of RNA-Seq data

# Installing packages from CRAN repo
install.packages('tidyverse')
install.packages('reshape2')
install.packages('dendextend')

# Installing packages from bioconductor repo
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("limma")
BiocManager::install("Biobase")

# Load libraries
library(edgeR)
library(DESeq2)
library(limma)
library(Biobase)
library(reshape2)
library(dendextend)


# LOADING DATA
# Import the meta data from the paper
SraRunTable.txt <- read.csv(file.choose())
SraRunTable <- SraRunTable.txt
SraRunTable

# Cleaning the column names
colnames(SraRunTable) = gsub("\\..exp.", " ", colnames(SraRunTable))

# Swapping out "-" and "+" for "minus" and "plus" because it will throw errors otherwise
SraRunTable$Cell_type = gsub("-", "_minus", SraRunTable$Cell_type)
# The "+" is a special character and should be escaped with "\\"
SraRunTable$Cell_type = gsub("\\+", "_plus", SraRunTable$Cell_type)
head(SraRunTable)

# Next we want to read in the data. Each sample counts are stored  in a seperate file into a count matrix
# Getting the list of all count files
file_list <- list.files(path = "../project/counts", full.names = T)

# Importation of files
library(readr)
GSM3690851_G01_htseq <- read_csv("counts/GSM3690851_G01_htseq.csv", 
                                 col_names = FALSE)

GSM3690852_G03_htseq <- read_csv("counts/GSM3690852_G03_htseq.csv", 
                                 col_names = FALSE)

GSM3690853_NG01_htseq <- read_csv("counts/GSM3690853_NG01_htseq.csv", 
                                  col_names = FALSE)

GSM3690854_NG03_htseq <- read_csv("counts/GSM3690854_NG03_htseq.csv", 
                                  col_names = FALSE)

GSM3690855_G05_htseq <- read_csv("counts/GSM3690855_G05_htseq.csv", 
                                 col_names = FALSE)

GSM3690856_G06_htseq <- read_csv("counts/GSM3690856_G06_htseq.csv", 
                                 col_names = FALSE)

GSM3690857_G07_htseq <- read_csv("counts/GSM3690857_G07_htseq.csv", 
                                 col_names = FALSE)

GSM3690858_G04_htseq <- read_csv("counts/GSM3690858_G04_htseq.csv", 
                                 col_names = FALSE)

GSM3690859_NG05_htseq <- read_csv("counts/GSM3690859_NG05_htseq.csv", 
                                  col_names = FALSE)

GSM3690860_NG06_htseq <- read_csv("counts/GSM3690860_NG06_htseq.csv", 
                                  col_names = FALSE)

GSM3690861_NG07_htseq <- read_csv("counts/GSM3690861_NG07_htseq.csv", 
                                  col_names = FALSE)

GSM3690862_NG04_htseq <- read_csv("counts/GSM3690862_NG04_htseq.csv", 
                                  col_names = FALSE)

GSM3690863_G10_htseq <- read_csv("counts/GSM3690863_G10_htseq.csv", 
                                 col_names = FALSE)

GSM3690864_G11_htseq <- read_csv("counts/GSM3690864_G11_htseq.csv", 
                                 col_names = FALSE)

GSM3690865_G12_htseq <- read_csv("counts/GSM3690865_G12_htseq.csv", 
                                 col_names = FALSE)

GSM3690866_G13_htseq <- read_csv("counts/GSM3690866_G13_htseq.csv", 
                                 col_names = FALSE)

GSM3690867_G14_htseq <- read_csv("counts/GSM3690867_G14_htseq.csv", 
                                 col_names = FALSE)

GSM3690868_G09_htseq <- read_csv("counts/GSM3690868_G09_htseq.csv", 
                                 col_names = FALSE)

GSM3690869_NG10_htseq <- read_csv("counts/GSM3690869_NG10_htseq.csv", 
                                  col_names = FALSE)

GSM3690870_NG11_htseq <- read_csv("counts/GSM3690870_NG11_htseq.csv", 
                                  col_names = FALSE)

GSM3690865_G12_htseq <- read_csv("counts/GSM3690865_G12_htseq.csv", 
                                 col_names = FALSE)

GSM3690866_G13_htseq <- read_csv("counts/GSM3690866_G13_htseq.csv", 
                                 col_names = FALSE)

GSM3690873_NG14_htseq <- read_csv("counts/GSM3690873_NG14_htseq.csv", 
                                  col_names = FALSE)

GSM3690874_NG09_htseq <- read_csv("counts/GSM3690874_NG09_htseq.csv", 
                                  col_names = FALSE)

GSM3690875_HC_G15_htseq <- read_csv("counts/GSM3690875_HC_G15_htseq.csv", 
                                    col_names = FALSE)

GSM3690876_HC_G16_htseq <- read_csv("counts/GSM3690876_HC_G16_htseq.csv", 
                                    col_names = FALSE)

GSM3690877_HC_G17_htseq <- read_csv("counts/GSM3690877_HC_G17_htseq.csv", 
                                    col_names = FALSE)

GSM3690878_HC_NG15_htseq <- read_csv("counts/GSM3690878_HC_NG15_htseq.csv", 
                                     col_names = FALSE)

GSM3690879_HC_NG16_htseq <- read_csv("counts/GSM3690879_HC_NG16_htseq.csv", 
                                     col_names = FALSE)

GSM3690880_HC_NG17_htseq <- read_csv("counts/GSM3690880_HC_NG17_htseq.csv", 
                                     col_names = FALSE)

GSM3690881_NS_G18_htseq <- read_csv("counts/GSM3690881_NS_G18_htseq.csv", 
                                    col_names = FALSE)

GSM3690882_NS_G19_htseq <- read_csv("counts/GSM3690882_NS_G19_htseq.csv", 
                                    col_names = FALSE)

GSM3690883_NS_G20_htseq <- read_csv("counts/GSM3690883_NS_G20_htseq.csv", 
                                    col_names = FALSE)

GSM3690884_NS_G21_htseq <- read_csv("counts/GSM3690884_NS_G21_htseq.csv", 
                                    col_names = FALSE)

GSM3690885_NS_NG18_htseq <- read_csv("counts/GSM3690885_NS_NG18_htseq.csv", 
                                     col_names = FALSE)

GSM3690886_NS_NG19_htseq <- read_csv("counts/GSM3690886_NS_NG19_htseq.csv", 
                                     col_names = FALSE)

GSM3690887_NS_NG20_htseq <- read_csv("counts/GSM3690887_NS_NG20_htseq.csv", 
                                     col_names = FALSE)

GSM3690888_NS_NG21_htseq <- read_csv("counts/GSM3690888_NS_NG21_htseq.csv", 
                                     col_names = FALSE)

counts = data.frame(GSM3690851_G01_htseq, GSM3690852_G03_htseq, GSM3690853_NG01_htseq, 
                    GSM3690854_NG03_htseq, GSM3690855_G05_htseq, GSM3690856_G06_htseq, 
                    GSM3690857_G07_htseq, GSM3690858_G04_htseq, GSM3690859_NG05_htseq, 
                    GSM3690860_NG06_htseq, GSM3690861_NG07_htseq, GSM3690862_NG04_htseq, 
                    GSM3690863_G10_htseq, GSM3690864_G11_htseq, GSM3690865_G12_htseq, 
                    GSM3690866_G13_htseq, GSM3690867_G14_htseq, GSM3690868_G09_htseq, 
                    GSM3690869_NG10_htseq, GSM3690870_NG11_htseq, GSM3690865_G12_htseq, 
                    GSM3690866_G13_htseq,GSM3690873_NG14_htseq, GSM3690874_NG09_htseq, 
                    GSM3690875_HC_G15_htseq, GSM3690876_HC_G16_htseq, GSM3690877_HC_G17_htseq, 
                    GSM3690878_HC_NG15_htseq, GSM3690879_HC_NG16_htseq, GSM3690880_HC_NG17_htseq, 
                    GSM3690881_NS_G18_htseq, GSM3690882_NS_G19_htseq, GSM3690883_NS_G20_htseq, 
                    GSM3690884_NS_G21_htseq, GSM3690885_NS_NG18_htseq, GSM3690886_NS_NG19_htseq, 
                    GSM3690887_NS_NG20_htseq, GSM3690888_NS_NG21_htseq)

# Extracting the GEO accession number for experiment identifier 
accession = gsub('^,*../project/counts/\\s*|\\_*', '', file_list)

# Reading in the gene list from the first count file
genes <- read.table(file_list[1], header = FALSE, sep = "\t") [,1]

# Reading in the counts from all the files
counts <- do.call(cbind, lapply(file_list, function(fn)read.table(fn, header = FALSE, sep = "\t" [, 2]))
colnames(counts) = accession
counts = data.frame(SYMBOL = genes, counts)

# Filter out the htseq stats
counts = counts[!c(grepl("_no_feature", counts$SYMBOL)|
                   grepl("_ambiguous", counts$SYMBOL)|
                  grepl("_too_low_aqual", counts$SYMBOL)|
                  grepl("_not_aligned", counts$SYMBOL)|
                  grepl("_alignment_not_unique", counts$SYMBOL)),]
tail(counts)
