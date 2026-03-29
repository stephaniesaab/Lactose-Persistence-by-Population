#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries
library(vcfR)
library(dplyr)
library(ggplot2)

#Need to merge VCF (Variant call format) data with metadata, want to align sampleIDs

#Read in VCF file
#Should output processed variant: 6776
vcf <- read.vcfR("../data/lactose_persistence.vcf", verbose = TRUE)

#Read metadata
metadata <- read.table("../data/20130606_g1k_3202_samples_ped_population.txt", header = TRUE, sep = " ") #Columns separated by spaces
names(metadata)

#Extract genotypes and sample IDs
#Extract genotype matrix (Converts genotypes into number of alternate alleles)
genotypes <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
nrow(genotypes) #6776, each row is a variant

#Transpose so samples are rows and SNPs are columns
genotypes_t <- t(genotypes)
ncol(genotypes_t) #6776

#Convert data to a dataframe to make sampleIDs a column
geno_df <- as.data.frame(genotypes_t)
geno_df$SampleID <- rownames(geno_df)

#Merge metadata with genotypes
geno_data <- geno_df %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))
nrow(geno_data) #3202 (3202 samples in 30X 1000 Genomes project)
ncol(geno_data) #6783 (6776 SampleIDs + SampleID + Sex + FamilyID + Population + Superpopulation + MotherID + FatherID )

#Check populations
print(table(geno_data$Population)) #Three letter codes for locations
print(table(geno_data$Superpopulation)) #Africa, America, East Asia, Europe, South Asia

