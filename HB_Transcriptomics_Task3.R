# Differential Expression Analysis of RNA-Seq data

# Installing packages from CRAN repo
install.packages('dplyr')
install.packages('tidyverse')

# Installing packages from bioconductor repo
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

# Load libraries
library(dplyr)
library(tidyverse)
library(readr)

# Importation of files
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
View(counts)

# Delete gene names that existed severally
dat = select(counts, -c(X1.1, X1.2, X1.3, X1.4, X1.5, X1.6, X1.7, X1.8, X1.9, X1.10,
                        X1.11, X1.12, X1.13, X1.14, X1.15, X1.16, X1.17, X1.18, X1.19,
                        X1.20, X1.21, X1.22, X1.23, X1.24, X1.25, X1.26, X1.27, X1.28,
                        X1.29, X1.30, X1.31, X1.32, X1.33, X1.34, X1.35, X1.36, X1.37))

View(dat)

# Get meta data
SraRunTable.txt <- read.csv(file.choose())
SraRunTable <- SraRunTable.txt
metadata <- SraRunTable

# Cleaning the column names
colnames(metadata) = gsub("\\..exp.", " ", colnames(metadata))

# Swapping out "-" and "+" for "minus" and "plus" because it will throw errors otherwise
metadata$Cell_type = gsub("-", "_minus", metadata$Cell_type)
# The "+" is a special character and should be escaped with "\\"
metadata$Cell_type = gsub("\\+", "_plus", metadata$Cell_type)
head(metadata)

# Get the columns needed. NOTE that column 2 will only be needed in the future for merging.
metadata.subset <- select(metadata, c(2,8,15,20,30))

# Swapping out "-" and "+" for "minus" and "plus" because it will throw errors otherwise
metadata.subset$Cell_type = gsub("-", "_minus", metadata.subset$Cell_type)

# The "+" is a special character and should be escaped with "\\"
metadata.subset$Cell_type = gsub("\\+", "_plus", metadata.subset$Cell_type)
head(metadata.subset)

# Lets perform addition operation on the datasets
metadata.modified <- metadata %>%
select(2,8,15,20,30)

# Change column name to sample and Substitute "RNA-Seq" with "X2"
colnames(metadata.modified) = gsub("Assay.Type", "samples", colnames(metadata.modified))
metadata.modified$samples = gsub("RNA-Seq", "X2", metadata.modified$samples)
head(metadata.modified)

# Lets take a look at what our gene expression data looks like
head(dat)

# reshaping data we will exclude the gene column. We are calling our sample column "sample" and their values "FPKM"
dat.long = dat %>%
  dplyr::rename(gene = X1) %>%
  gather(key = 'samples', value = 'FPKM', -gene)

# Filter out the htseq stats
dat.long = dat.long[!c(grepl("_no_feature", dat.long$gene)|
                     grepl("_ambiguous", dat.long$gene)|
                     grepl("_too_low_aqual", dat.long$gene)|
                     grepl("_not_aligned", dat.long$gene)|
                     grepl("_alignment_not_unique", dat.long$gene)),]

# Let us replace all NA values with 0
dat.long[is.na(dat.long)] <- 0

# join dataframes = dat.long + metadata.modified
dat.long = dat.long %>%
  left_join(., metadata.modified, by = c("samples"))
View(dat.long)


#######################VISUALIZATION#####################################

# loading the library
library(ggplot2)
# 1. Bar plot
dat.long %>%
  filter(FPKM >= 50) %>%
  filter(!is.na(Treatment)) %>%
  ggplot(., aes(x = Treatment, y = FPKM, fill = Mouse_ID)) +
  geom_col()

# 2. Density plot
dat.long %>%
  filter(FPKM >= 50) %>%
  filter(!is.na(Treatment)) %>%
  ggplot(., aes(x = Cell_type, fill = Treatment)) +
  geom_density(alpha = 0.3)

# 3. Heat map
mouse.of.interest <- c("MouseNr 01","MouseNr 02","MouseNr 03","MouseNr 04",
                       "MouseNr 05","MouseNr 06","MouseNr 07","MouseNr 08",
                       "MouseNr 09","MouseNr 10","MouseNr 11","MouseNr 12",
                       "MouseNr 13","MouseNr 14","MouseNr 15","MouseNr 16",
                       "MouseNr 17","MouseNr 18","MouseNr 19","MouseNr 20",
                       "MouseNr 21")
dat.long %>%
  filter(Mouse_ID%in%mouse.of.interest) %>%
  filter(FPKM >= 50) %>%
  filter(!is.na(Mouse_ID)) %>%
  ggplot(., aes(x = Treatment, y = Mouse_ID, fill = FPKM)) +
  geom_tile()

