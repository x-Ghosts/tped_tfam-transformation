install.packages("tidyverse")

library(tidyverse)
dataset <- read.table("Desktop/dataset_editing/angora/angora_plos.txt", header = FALSE, stringsAsFactors = FALSE)
View(dataset)
colnames(dataset) <- c("FID", "IID", "SNP", "A1", "A2")

#manifest
manifest <- readLines("Desktop/dataset_editing/manifest_illumina_2.csv")
headers <- grep("^\\[Assay\\]", manifest)
manifest <- read.csv("Desktop/dataset_editing/manifest_illumina_2.csv", skip = headers, stringsAsFactors = FALSE)
View(manifest)


# Making tfam - FID IID PID MID SEX PHENO

tfam <- dataset %>%
  select(FID,IID) %>%
  distinct() %>%
  mutate(PID=0, MID=0, SEX=0, PHENO=-9)

View(tfam)

write.table(tfam,file = "Desktop/dataset_editing/angora/angora_plos.tfam", row.names = FALSE, col.names = FALSE, quote = FALSE)
View(dataset);View(tfam)


# Making TPED - 

## Format:      tped
## Description: Output to TPED format with the first four columns chr name gen_pos
## pos, and the rest for genotypes. Variant tools cannot import from
## this format because it does not contain information about reference
## genome.

## Columns:
## 1            chromosome (without leading chr)
## 2             Locus name
## 3            Genetic distance, left empty
## 4            Physical position
## 5            genotype

## First of all, start a dataframe with all chromosomes from the manifest

#tped_file <- data.frame(manifest$Chr)
#View(tped_file)

matched_snps <- manifest %>% filter(manifest$Name %in% dataset$SNP)
View(matched_snps)

tped_file <- data.frame(matched_snps$Name)

matched_chr <- manifest %>% 
  filter(Name %in% tped_file$matched_snps.Name) %>%
  select(Name, Chr)

matched_chr <- data.frame(matched_chr$Chr, matched_chr$Name)
tped_file <- matched_chr


matched_physical_position <- manifest %>% 
  filter(Name %in% tped_file$matched_chr.Name) %>%
  select(Name, Chr, MapInfo)

tped_file <- matched_physical_position
tped_file <- tped_file %>%
  mutate(Genomic_Position=0)

tped_file <- data.frame(tped_file$Chr, tped_file$Name, tped_file$Genomic_Position, tped_file$MapInfo)
View(tped_file)




geno <- dataset %>%
  mutate(FID = paste(FID, IID, sep = "_"))
View(geno)

tped_file_wide <- dataset %>%
  select(SNP, IID, A1, A2) %>%
  pivot_wider(names_from = IID, values_from = c(A1, A2))

view(tped_file_wide)

samples <- unique(dataset$IID)
ordered_cols <- c("SNP", unlist(lapply(samples, function(id) c(paste0("A1_", id), paste0("A2_", id)))))
tped_file_wide <- tped_file_wide[, ordered_cols]


## This is the FINAL PART OF THE MERGE, TO BE CONTINUED BY FRIDAY MORNING 

#tped_file <- tped_file_wide %>%
#  filter(SNP %in% tped_file$Name) %>%
#  select(tped_file$tped_file.Chr, tped_file$tped_file.Name, tped_file$tped_file.Genomic_Position, tped_file$tped_file.MapInfo, ordered_cols)

tped_file_full <- tped_file_wide %>%
  filter(SNP %in% tped_file$tped_file.Name) %>%
  left_join(tped_file, by = c("SNP" = "tped_file.Name")) %>%
  select(Chr, SNP, Genomic_Position, MapInfo, all_of(ordered_cols))
#####

# First Step
#matched_snps <- manifest %>% filter(Name %in% dataset$SNP)

# Step 2: Keep only relevant columns for TPED
#tped_file <- matched_snps %>%
#  select(Chr, Name, MapInfo) %>%
#  mutate(Genomic_Position = 0) %>%
#  select(Chr, Name, Genomic_Position, MapInfo)

View(tped_file)


