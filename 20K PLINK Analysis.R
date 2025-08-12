#This script integrates 20K SNP data from Howard et al. 2021, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#PLINK files for JD-PFR analysis previously are reused, and the 20K data formatted as close to PLINK format as possible in Excel first

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs")

#Run PLINK to extract SNPs and recode alleles
system("plink --file JKI_PLINK --extract 50K_extract --recode A --out 20K_JD_PFR")