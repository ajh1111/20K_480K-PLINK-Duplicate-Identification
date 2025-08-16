#This script integrates 20K SNP data from Larsen et al. 2025, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#PLINK files for JD-PFR analysis previously are reused.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Convert Larsen et al. 2025 data from VCF to PLINK .ped and .map
system("plink --vcf Danish_20K.vcf --const-fid 0 --allow-extra-chr --recode tab --out Danish_20K")

#Use PLINK to extract overlapping SNPs from JD and PFR samples
system("plink --file JD_PFR_All --extract JD_PFR_Danish_Dutch_ExtractList.txt --make-bed --out Danish_Dutch_JD_PFR")
system("plink --bfile Danish_Dutch_JD_PFR --recode --tab --out Danish_Dutch_JD_PFR")

#Use PLINK to extract overlapping SNPs from 20K Larsen et al. 2025 samples, recode SNP positions, and output as .ped and .map
system("plink --file Danish_20K --extract Danish_Dutch_ExtractList.txt --allow-extra-chr --update-chr Danish_Dutch_chr.txt --update-cm Danish_Dutch_cm.txt --update-map Danish_Dutch_map.txt --make-bed --out Danish_20K")
system("plink --bfile Danish_20K  --update-name Danish_Dutch_name.txt --recode tab --out Danish_Ready")