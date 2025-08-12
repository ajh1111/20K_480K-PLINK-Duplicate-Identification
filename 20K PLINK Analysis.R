#This script integrates 20K SNP data from Howard et al. 2021, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#PLINK files for JD-PFR analysis previously are reused, and the 20K data formatted as close to PLINK format as possible in Excel first

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Initial JD_PFR .ped file curation

#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Raw.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Remove triploids from .ped file
ids_to_remove <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_remove <- ids_to_remove$V1
ped <- ped[!(ped$V1 %in% ids_to_remove), ]

#Add empty columns for PLINK formatting
ped <- add_column(ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
ped <- add_column(ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_All.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Initial JD_PFR .map file curation

#Load .map file and BLAST results
map <- read.csv("JD_PFR_Raw.map", header = FALSE, sep ="\t")
BLAST <- read.csv("BLAST results.tsv", header = FALSE, sep = "\t")

#Match Marker IDs
match_ids <- match(map$V2, BLAST$V2)
matched <- !is.na(match_ids)
map[matched, ] <- BLAST[match_ids[matched], ]

#Remove header row, and set any chromosome numbers larger than 17 to 0
map = map[-1, ]
map$V1 <- as.numeric(as.character(map$V1))
map[map$V1 > 17, c("V1", "V4")] <- 0

#Save .map file (with locations)
write.table(map, "JD_PFR_All.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Run PLINK to extract SNPs and recode alleles
system("plink --file JD_PFR_All --update-alleles Recode_50K.txt  --extract 50K_extract.txt --make-bed --out 20K_JD_PFR")
system("plink --bfile 20K_JD_PFR --recode --tab --out 20K_JD_PFR")


##Combining Howard et al. 2021 20K data with JD_PFR samples

#clear workspace
rm(list=ls())

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Load Howard 20K data and transpose
HW <- read.delim("Howard_2021_20K.txt", header = FALSE, stringsAsFactors = FALSE)
HW_t <- as.data.frame(t(HW))

#Add empty columns for PLINK formatting
HW_t <- add_column(HW_t, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
HW_t <- add_column(HW_t, Fid = 0, .before = "V1" )

#Remove header row 
HW_t = HW_t[-1, ]

#Remove row and column names
colnames(HW_t) <- NULL
rownames(HW_t) <- NULL

#Load in PLINK .ped
ped <- read.csv("20K_JD_PFR.ped", header = FALSE,sep = "\t")

#Bind JKI data onto PLINK .ped
names(HW_t) <- names(ped)
combined_ped <- bind_rows(ped, HW_t)

#Load .map file
map <- read.csv("20K_JD_PFR.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "20K_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "20K_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
