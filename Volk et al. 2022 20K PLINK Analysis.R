#This script integrates 20K SNP data from Volk et al. 2022, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#PLINK files for JD-PFR analysis previously are reused. SNP list for Larsen et al. data is reused too.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")


##Convert files, format, extract chosen SNPs.

#Convert Volk et al. 2022 data from VCF to PLINK .ped and .map
system("plink --vcf Volk_20K.vcf --const-fid 0 --allow-extra-chr --recode tab --out Volk_20K")

#Use PLINK to extract overlapping SNPs from JD and PFR samples (using list from Larsen et al. R script)
system("plink --file JD_PFR_All --extract JD_PFR_Danish_Dutch_ExtractList.txt --make-bed --out Volk_JD_PFR")
system("plink --bfile Volk_JD_PFR --recode --tab --out Volk_JD_PFR")

#Use PLINK to extract overlapping SNPs from 20K Larsen et al. 2025 samples, recode SNP positions, and output as .ped and .map
system("plink --file Volk_20K --extract Danish_Dutch_ExtractList.txt --allow-extra-chr --update-chr Danish_Dutch_chr.txt --update-cm Danish_Dutch_cm.txt --update-map Danish_Dutch_map.txt --make-bed --out Volk_20K")
system("plink --bfile Volk_20K  --update-name Danish_Dutch_name.txt --recode tab --out Volk_Ready")


##Combining Volk et al. 2022 20K data with JD_PFR samples

#clear workspace
rm(list=ls())

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Load JD_PFR data
JD_ped <- read.csv("Volk_JD_PFR.ped", header = FALSE,sep = "\t")

#Load Volk et al. 2022 20K data, add prefix to names, set irrelevant tags to 0
VK_ped <- read.csv("Volk_Ready.ped", header = FALSE,sep = "\t")
VK_ped$V2 <- paste0("VK_", VK_ped$V2) 
VK_ped$V1 <- 0
VK_ped$V3 <- 0
VK_ped$V4 <- 0

#Bind Volk et al. 2022 20K data onto PLINK .ped, set other variables to 0
names(VK_ped) <- names(JD_ped)
combined_ped <- bind_rows(JD_ped, VK_ped)

#Load .map file
map <- read.csv("Volk_JD_PFR.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "Volk_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "Volk_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##Running PLINK Duplicate Analysis

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Run PLINK
system("plink --file Volk_PLINK --missing-genotype 0 --genome full ")

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/Volk et al. 2022 20K Results/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##Grouping duplicates

#Filter for PI_HAT >0.96 (duplicate threshold)
genome <- genome[!(genome$PI_HAT < 0.96), ]
genome <- subset(genome, select = c("IID1","IID2"))

#Group duplicates with igraph
graph <- graph_from_data_frame(genome, directed = FALSE)
components <- components(graph)

#Sort groupings by number of duplicates
group_sizes <- table(components$membership)
sorted_group_ids <- order(group_sizes)
new_ids <- match(components$membership, sorted_group_ids)
V(graph)$group <- new_ids
grouped_samples <- split(names(components$membership), new_ids)

#Pad group with length less than max length with NA's
max_len <- max(sapply(grouped_samples, length))
padded_list <- lapply(grouped_samples, function(x) {c(x, rep(" ", max_len - length(x)))})

#Write groupings to a dataframe
dd <- as.data.frame(do.call(rbind, padded_list))

#Add a number for each group
dd <- cbind(Group = seq_len(nrow(dd)), dd)

# Add the number of duplicates in each grouping
sample_counts <- rowSums(dd[, -1] != " ")
dd <- add_column(dd, SampleCount = sample_counts, .after = "Group")

#Rename columns
colnames(dd) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13")

#Save .csv of duplicate groupings
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/Volk et al. 2022 20K Results/Grouped_Duplicates.csv", row.names = FALSE)

