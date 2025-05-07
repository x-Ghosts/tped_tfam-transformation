# library to import : 

library(tidyverse)
library(tidyr)



dataset <- read.table("angora_plos.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(dataset) <- c("FID", "IID", "SNP", "A1", "A2")


manifest <- read.csv("manifest_illumina_2.csv", sep = ",", skip = 8, header = TRUE, stringsAsFactors = FALSE)
header_manifest <- read.csv("manifest_illumina_2.csv", skip = 7, header = FALSE, nrows = 1)
colnames(manifest) <- header_manifest


# MAP build :

list_SNP <- manifest %>% filter(Name %in% dataset$SNP)
map_file <- list_SNP %>%
  mutate(Morgan_distance = 0 ) %>%
  select(Chr,Name,Morgan_distance,MapInfo)
write.table("angora_plos.map", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# PED build
subset_dataset <- subset(dataset, dataset$IID == "1_1142_A")
list_of_ordered_snp <- map_file$Name

ped_file <- subset_dataset %>%
  mutate(PID=0, MID=0, sex=0, pheno=-9) %>%
  pivot_wider(id_cols = c(FID,IID,PID,MID,sex,pheno), names_from = SNP, values_from = c(A1,A2))

ped_ordered_cols <- unlist(lapply(subset_dataset$SNP, function(SNP) c(paste0("A1_", SNP), paste0("A2_", SNP))))
ped_file <- ped_file %>% select(FID,IID,PID,MID,sex,pheno, all_of(ped_ordered_cols))

colnames(subset_dataset)
