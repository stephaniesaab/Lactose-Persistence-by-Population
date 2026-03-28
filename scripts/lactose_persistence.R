#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries
library(vcfR)
library(dplyr)
library(ggplot2)

#Need to merge VCF (Variant call format) data with metadata, want to align sampleIDs

#Extract genotypes into a matrix and then join matrix by labels

#Read in VCF file
vcf <- read.vcfR("../data/lactose_data.vcf")

#Read metadata
metadata <- read.table("../data/20130606_g1k_3202_samples_ped_population.txt", header = TRUE, sep = "\t")
