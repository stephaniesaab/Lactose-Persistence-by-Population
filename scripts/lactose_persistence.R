#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries ====
library(vcfR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(randomForest)
library(pegas)

#TODO =====
# 1) Fix the EDA (take out outliers part, i dont think that makes sense for SNP data)
# 2) Random forest model
# 3) K means model
# 4) Figures for PCA
# 5) Figures for K means
# 6) Figures for RF
#Need to merge VCF (Variant call format) data with metadata, want to align sampleIDs

#Read in VCF file ====
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
geno_df <- as.data.frame(genotypes_t) #nrow  = 
geno_df$SampleID <- rownames(geno_df)
geno_df <- geno_df %>%
  filter(SampleID %in% metadata$SampleID)

#Merge metadata with genotypes
geno_data <- geno_df %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))
nrow(geno_data) #3202 (3202 samples in 30X 1000 Genomes project)
ncol(geno_data) #6783 (6776 SampleIDs + SampleID + Sex + FamilyID + Population + Superpopulation + MotherID + FatherID )

#Check populations
print(table(geno_data$Population)) #Three letter codes for locations
print(table(geno_data$Superpopulation)) #Africa, America, East Asia, Europe, South Asia


#Write out CSV with variants and data
#write.csv(geno_data, "../data/genotype_data.csv")

# EDA ====
#Exploratory data analysis to identify which subpopulations to extract before doing clustering by genotypes

##1. Initial Subsetting and cleaning ====

#Filter for target populations first (AFR, EUR, EAS, SAS).
#Ensures HWE and MAF calcs are for specific analysis groups.

geno_subset <- geno_data %>% 
  filter(Superpopulation %in% c("AFR", "EUR", "EAS", "SAS"))

#Separate genotype columns from the metadata columns
metadata_cols <- c("SampleID", "FamilyID", "MotherID", "FatherID", "Sex", "Population", "Superpopulation")
all_snps <- setdiff(names(geno_subset), metadata_cols)
length(all_snps) #= 6776 SNPs

##2. SNP Filtering (MAF and informative) ====
#Remove rare variants and monomorphic loci
#Calculate the sum of alleles for each SNP
# A sum of 0 means the variant (minor allele) is not present in the 3,202 samples
mat_temp <- as.matrix(geno_subset[, all_snps])
snp_sums <- colSums(mat_temp)
n_samples <- nrow(geno_subset) #2712
snp_sds <- apply(mat_temp, 2, sd, na.rm = TRUE)

#How many values are missing:
any(is.na(mat_temp)) #FALSE so no missing values

#Remove SNPs where everyone has the variant (sum == 2 * number of samples (3202))
#Remove SNPs that are rare in the population (<0.001)
keep_snps <-  snp_sums < (2*n_samples*(1- 0.001)) &
  snp_sums > (0.001*n_samples) &
  snp_sds > 0
active_snps <- names(which(keep_snps))
keep_snps[is.na(keep_snps)] <- FALSE

#Subset data to keep only informative SNPs
geno_final <- geno_subset[, 
                          c("SampleID", "Sex", "Population", "Superpopulation",
                            active_snps)]
length(active_snps) #2516
geno_matrix <- as.matrix(geno_final[, active_snps])
rownames(geno_matrix) <- geno_final$SampleID
head(geno_matrix[, 1:5])

print(paste("Final SNPs for analysis:", length(active_snps))) # Should be 2556

##3. HWE filtering ====
# 1 = p^2 + 2pq + q^2
#Identify SNPs that deviate from HWE within specific groups
#Prevents technical errors from biasing PCA and Random forest
filter_hwe_by_pop <- function(df, pop_col = "Superpopulation", p_threshold = 0.05) {
  
  #Get columns that are SNPs
  snp_cols <- setdiff(names(df), c("SampleID", "FamilyID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"))
  
  #Split dataframe by Superpopulation labels
  super_pops <- split(df, df[[pop_col]])
  
  #Create a list to track SNPs that fail in any population
  snps_to_rm <- c()
  
  for (pop_name in names(super_pops)) {
    message(paste("Checking HWE for population:", pop_name))
    
    #Get the genotype matrix for this pop
    sub_matrix <- super_pops[[pop_name]][, snp_cols]
    
    #Chi square test (Pearson's X^2 test)
    for (snp in snp_cols) {
      genotypes <- sub_matrix[[snp]]
      
      #Get genotype counts: 0 (Homoz-Ref), 1(Het), 2(Homoz-Alt)
      n0 <- sum(genotypes == 0)
      n1 <- sum(genotypes == 1)
      n2 <- sum(genotypes ==2)
      n <- n0 + n1 + n2
      
      if (n > 0) {
        #Calculate allele frequencies
        p <- (2*n0 + n1) / (2*n)
        q <- 1 - p
        
        #Expected counts in HWE
        exp0 <- (p^2) * n
        exp1 <- (2*p*q) * n
        exp2 <- (q^2) * n
        
        #Chi square test
        #Avoid division by zero if expectations are very small
        if (all(c(exp0, exp1, exp2) > 5)) {
          chisq <- ((n0 - exp0)^2 / exp0) + ((n1 - exp1)^2 / exp1) + ((n2 - exp2)^2 / exp2)
          p_val <- pchisq(chisq, df = 1, lower.tail = FALSE)
          
          if (p_val < p_threshold) {
            snps_to_rm <- c(snps_to_rm, snp)
          }
        }
      }
    }
  }
  
  #Return unique list of SNPs that failed HWE in at least one group
  return(unique(snps_to_rm))
}

# Identify SNPs that fail HWE
# hwe_fails <- filter_hwe_by_pop(geno_data_filt, pop_col = "Superpopulation")
# length(hwe_fails) #Remove 679 SNPs

# #Remove the SNPs
# geno_subset_filt <- geno_data_filt %>% 
#   select(-all_of(hwe_fails))
# snp_cols_filt <- genotypes_filt %>% 
#   select(-all_of(hwe_fails))
# ncol(snp_cols_filt) #2833

##3. QC ====
# Population Sample Size Bar Chart
pop_counts <- geno_final%>%
  count(Superpopulation, Population)

ggplot(pop_counts, aes(x = reorder(Population, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Sample Size by Population",
    x = "Population",
    y = "Sample Count",
    fill = "Superpopulation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Population Counts by Superpopulation
ggplot(pop_counts %>% group_by(Superpopulation) %>% summarise(n = sum(n)),
       aes(x = reorder(Superpopulation, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  labs(title = "Sample Size by Superpopulation", x = "Superpopulation", y = "Count") +
  theme_minimal()

## Per-Sample Heterozygosity ====
het_rates <- rowMeans(geno_final[, active_snps] == 1)
het_df <- data.frame(
  SampleID        = geno_subset$SampleID,
  Heterozygosity  = het_rates,
  Superpopulation = geno_subset$Superpopulation
)

het_mean <- mean(het_df$Heterozygosity)
het_sd   <- sd(het_df$Heterozygosity)
het_df$Outlier <- abs(het_df$Heterozygosity - het_mean) > 3 * het_sd

cat("Heterozygosity outliers:", sum(het_df$Outlier), "\n")
print(het_df[het_df$Outlier, ])

ggplot(het_df, aes(x = Superpopulation, y = Heterozygosity, fill = Superpopulation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Outlier), width = 0.2, size = 0.8, alpha = 0.4) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
  geom_hline(yintercept = c(het_mean - 3*het_sd, het_mean + 3*het_sd),
             linetype = "dashed", color = "red", linewidth = 0.6) +
  labs(
    title = "Per-Sample Heterozygosity by Superpopulation",
    subtitle = "Red dashed lines = mean ± 3 SD thresholds",
    x = "Superpopulation",
    y = "Heterozygosity Rate"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


#PCA =====
#Prepare numeric genotypes (Isolate only the SNP columns for calcs)

#Normalization
#Using scale. = TRUE in prcomp to standardize data
#Scales each SNP to have unit variance (mean = 0, SD = 1)

#Perform PCA
pca_result <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)
#Summary results of PCA
summary_pca <- summary(pca_result)

#Variance explanation (Scree plot)
pca_var <- pca_result$sdev^2
var_explained <- pca_var / sum(pca_var)
cumu_var <- cumsum(var_explained)

#Figures for PCA: Scree plot
plot(var_explained[1:10], type = "b",
     main = "Figure 1: Scree Plot of LCT and MCM6 Region PCA",
     xlab = "Principal Component",
     ylab = "Proportion of Variance",
     col = "blue", pch = 19)

pca_df <- as.data.frame(pca_result$x)
pca_df$Superpopulation <- geno_final$Superpopulation

pca_var     <- pca_result$sdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 1)

pca_df <- data.frame(
  PC1             = pca_result$x[, 1],
  PC2             = pca_result$x[, 2],
  PC3             = pca_result$x[, 3],
  Superpopulation = geno_final$Superpopulation,
  Population      = geno_final$Population
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA of Genotype Data (pre-HWE filtering)",
    subtitle = "Colored by Superpopulation",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC2 (", pca_var_pct[2], "% variance)")
  ) +
  theme_minimal()

ggplot(pca_df, aes(x = PC1, y = PC3, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA: PC1 vs PC3 (pre-HWE filtering)",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC3 (", pca_var_pct[3], "% variance)")
  ) +
  theme_minimal()

scree_df <- data.frame(PC = 1:20, VariancePct = pca_var_pct[1:20])
ggplot(scree_df, aes(x = PC, y = VariancePct)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot (pre-HWE filtering)", x = "Principal Component", y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:20) +
  theme_minimal()

#Identifying top SNPs ====
#Goal is to identify lactose persistence across genetic groups so want to colour PCA by superpopulation

#PCA of lactose persistence region
ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal()+
  labs(title = "Figure 2: PCA of Lactose Persistence by Region",
       subtitle = "Samples clustered by superpopulation",
       x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"))


#Identify top contributing SNPs (Loadings) ====
snp_loadings <- data.frame(
  SNP = rownames(pca_result$rotation),
  PC1_loading = abs(pca_result$rotation[,1])
)

top_snps <- snp_loadings %>% 
  arrange(desc(PC1_loading)) %>% 
  head(10)
#Get top 10 SNPs
print(top_snps)


# 5. Variant Density Across the Region ====
#Define slice coordinates
region_start <- 135700000
region_end   <- 136100000

#Filter the SNP positions to only include this slice
snp_positions_slice <- snp_positions_filt %>%
  filter(POS >= region_start & POS <= region_end)

# Plot the specific region
ggplot(snp_positions_slice, aes(x = POS)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  labs(
    title = "Figure 3: Variant Density Across LCT/MCM6 Region",
    subtitle = paste0("Genomic Slice: ", scales::comma(region_start), " - ", scales::comma(region_end), " bp"),
    x = "Genomic Position (Chromosome 2)",
    y = "Number of Variants"
  ) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal()

#K-means ====

#Get optimal clusters
#Expected 3 or 5 based on Superpopulations (AFR, EUR, EAS, SAS)
set.seed(1516)
#Run K-means only on the numeric SNP columns
#Use active_snps to ensure we aren't sending 'SampleID' to the math function
km_res <- kmeans(geno_matrix, centers = 3, nstart = 25)

#Add cluster assignments for plotting
geno_final$Cluster <- as.factor(km_res$cluster)

#Check results against the actual populations
table(geno_final$Superpopulation, geno_final$Cluster)

#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries ====
library(vcfR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(randomForest)
library(pegas)
#Need to merge VCF (Variant call format) data with metadata, want to align sampleIDs

setwd("~/BINF_DESKTOP/BINF6970/assignments/assignment03/Lactose-Persistence-by-Population/scripts")

#Read in VCF file ====
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
geno_df <- as.data.frame(genotypes_t) #nrow  = 
geno_df$SampleID <- rownames(geno_df)
geno_df <- geno_df %>%
  filter(SampleID %in% metadata$SampleID)

#Merge metadata with genotypes
geno_data <- geno_df %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))
nrow(geno_data) #3202 (3202 samples in 30X 1000 Genomes project)
ncol(geno_data) #6783 (6776 SampleIDs + SampleID + Sex + FamilyID + Population + Superpopulation + MotherID + FatherID )

#Check populations
print(table(geno_data$Population)) #Three letter codes for locations
print(table(geno_data$Superpopulation)) #Africa, America, East Asia, Europe, South Asia


#Write out CSV with variants and data
#write.csv(geno_data, "../data/genotype_data.csv")

# EDA ====
#Exploratory data analysis to identify which subpopulations to extract before doing clustering by genotypes

#Separate genotype columns from the metadata columns
genotype_columns <- geno_data %>% #ncol = 6776 variants
  select(-SampleID, -FamilyID, -FatherID, -MotherID, -Sex, -Population, -Superpopulation)

#How many values are missing:
missing_snps <- is.na(genotype_columns)
sum(missing_snps) #0 so no missing values

#Calculate the sum of alleles for each SNP
# A sum of 0 means the variant (minor allele) is not present in the 3,202 samples
snp_sums <- colSums(genotype_columns, na.rm = TRUE)

#Identify SNPs with at least one variant (sum > 0)
#Remove SNPs where everyone has the variant (sum == 2 * number of samples (3202))
#Remove SNPs that are rare in the population (<0.001)
keep_snps <- snp_sums > 0 & snp_sums < 2*nrow(geno_data) & snp_sums > 0.001*nrow(geno_data)

#Subset data to keep only informative SNPs
geno_data_filt <- geno_data[, c(names(which(keep_snps)), "SampleID", "Population", "Superpopulation", "Sex")]
genotypes_filt <- genotype_columns[, names(which(keep_snps))]

# ***************** EDA added by Maryanne********************************
#Making the SNP columns available
snp_cols <- geno_data_filt %>%
  select(-SampleID, -Sex, -Population, -Superpopulation) %>%
  names()

geno_matrix <- geno_data_filt %>%
  select(all_of(snp_cols)) %>%
  as.matrix()
rm(geno_matrix)

# 1. Population Sample Size Bar Chart ====
pop_counts <- geno_data_filt %>%
  count(Superpopulation, Population)

ggplot(pop_counts, aes(x = reorder(Population, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Sample Size by Population",
    x = "Population",
    y = "Sample Count",
    fill = "Superpopulation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Population Counts by Superpopulation
ggplot(pop_counts %>% group_by(Superpopulation) %>% summarise(n = sum(n)),
       aes(x = reorder(Superpopulation, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  labs(title = "Sample Size by Superpopulation", x = "Superpopulation", y = "Count") +
  theme_minimal()

# 2. Minor Allele Frequency (MAF) Distribution ====
maf_values <- sapply(snp_cols, function(snp) {
  allele_freq <- mean(geno_data_filt[[snp]], na.rm = TRUE) / 2
  return(min(allele_freq, 1 - allele_freq))
})

maf_df <- data.frame(SNP = names(maf_values), MAF = maf_values)

ggplot(maf_df, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, label = "MAF = 0.05", vjust = 2, color = "red", size = 3.5) +
  labs(
    title = "Minor Allele Frequency Distribution",
    subtitle = paste0("n = ", nrow(maf_df), " SNPs after MAF filtering, before HWE filtering"),
    x = "Minor Allele Frequency",
    y = "Number of SNPs"
  ) +
  theme_minimal()

#COME BACK HERE AND FIND OUT WHERE GENO M<ATRIX WENT 
## 3. Per-Sample Heterozygosity ====
het_rates <- apply(geno_matrix, 1, function(x) {
  sum(x == 1, na.rm = TRUE) / sum(!is.na(x))
})

het_df <- data.frame(
  SampleID        = geno_data_filt$SampleID,
  Heterozygosity  = het_rates,
  Superpopulation = geno_data_filt$Superpopulation
)

het_mean <- mean(het_df$Heterozygosity)
het_sd   <- sd(het_df$Heterozygosity)
het_df$Outlier <- abs(het_df$Heterozygosity - het_mean) > 3 * het_sd

cat("Heterozygosity outliers:", sum(het_df$Outlier), "\n")
print(het_df[het_df$Outlier, ])

ggplot(het_df, aes(x = Superpopulation, y = Heterozygosity, fill = Superpopulation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Outlier), width = 0.2, size = 0.8, alpha = 0.4) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
  geom_hline(yintercept = c(het_mean - 3*het_sd, het_mean + 3*het_sd),
             linetype = "dashed", color = "red", linewidth = 0.6) +
  labs(
    title = "Per-Sample Heterozygosity by Superpopulation",
    subtitle = "Red dashed lines = mean ± 3 SD thresholds",
    x = "Superpopulation",
    y = "Heterozygosity Rate"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 4. PCA Plot ====
# Impute any NAs with column means before running PCA
geno_matrix_imputed <- apply(geno_matrix, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

pca_result <- prcomp(geno_matrix_imputed, center = TRUE, scale. = FALSE)

pca_var     <- pca_result$sdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 1)

pca_df <- data.frame(
  PC1             = pca_result$x[, 1],
  PC2             = pca_result$x[, 2],
  PC3             = pca_result$x[, 3],
  Superpopulation = geno_data_filt$Superpopulation,
  Population      = geno_data_filt$Population
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA of Genotype Data (pre-HWE filtering)",
    subtitle = "Colored by Superpopulation",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC2 (", pca_var_pct[2], "% variance)")
  ) +
  theme_minimal()

ggplot(pca_df, aes(x = PC1, y = PC3, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA: PC1 vs PC3 (pre-HWE filtering)",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC3 (", pca_var_pct[3], "% variance)")
  ) +
  theme_minimal()

scree_df <- data.frame(PC = 1:20, VariancePct = pca_var_pct[1:20])
ggplot(scree_df, aes(x = PC, y = VariancePct)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot (pre-HWE filtering)", x = "Principal Component", y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:20) +
  theme_minimal()
# 5. Variant Density Across the Region ====
snp_positions <- data.frame(
  CHROM = getCHROM(vcf),
  POS   = as.numeric(getPOS(vcf)),
  ID    = getID(vcf)
)

snp_positions_filt <- snp_positions %>%
  filter(ID %in% snp_cols)

ggplot(snp_positions_filt, aes(x = POS)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  labs(
    title = "Variant Density Across the Lactose Persistence Locus",
    subtitle = "After MAF filtering, before HWE filtering",
    x = "Genomic Position (bp)",
    y = "Number of Variants"
  ) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal()

if (length(unique(snp_positions_filt$CHROM)) > 1) {
  ggplot(snp_positions_filt, aes(x = POS)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    facet_wrap(~CHROM, scales = "free_x") +
    labs(title = "Variant Density by Chromosome", x = "Position", y = "Variants") +
    scale_x_continuous(labels = scales::comma) +
    theme_minimal()
}
# **********************************************
#Filter for geographic locations that will likely show differences in lactose persistence based on literature review

# HWE Filtering Function by Superpopulation ====
# 1 = p^2 + 2pq + q^2
#Identify SNPs that deviate from HWE within specific groups
#Prevents technical errors from biasing PCA and Random forest
filter_hwe_by_pop <- function(df, pop_col = "Superpopulation", p_threshold = 0.05) {
  
  #Get columns that are SNPs
  snp_cols <- setdiff(names(df), c("SampleID", "FamilyID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"))
  
  #Split dataframe by Superpopulation labels
  super_pops <- split(df, df[[pop_col]])
  
  #Create a list to track SNPs that fail in any population
  snps_to_rm <- c()
  
  for (pop_name in names(super_pops)) {
    message(paste("Checking HWE for population:", pop_name))
    
    #Get the genotype matrix for this pop
    sub_matrix <- super_pops[[pop_name]][, snp_cols]
    
    #Chi square test (Pearson's X^2 test)
    for (snp in snp_cols) {
      genotypes <- sub_matrix[[snp]]
      
      #Get genotype counts: 0 (Homoz-Ref), 1(Het), 2(Homoz-Alt)
      n0 <- sum(genotypes == 0)
      n1 <- sum(genotypes == 1)
      n2 <- sum(genotypes ==2)
      n <- n0 + n1 + n2
      
      if (n > 0) {
        #Calculate allele frequencies
        p <- (2*n0 + n1) / (2*n)
        q <- 1 - p
        
        #Expected counts in HWE
        exp0 <- (p^2) * n
        exp1 <- (2*p*q) * n
        exp2 <- (q^2) * n
        
        #Chi square test
        #Avoid division by zero if expectations are very small
        if (all(c(exp0, exp1, exp2) > 5)) {
          chisq <- ((n0 - exp0)^2 / exp0) + ((n1 - exp1)^2 / exp1) + ((n2 - exp2)^2 / exp2)
          p_val <- pchisq(chisq, df = 1, lower.tail = FALSE)
          
          if (p_val < p_threshold) {
            snps_to_rm <- c(snps_to_rm, snp)
          }
        }
      }
    }
  }
  
  #Return unique list of SNPs that failed HWE in at least one group
  return(unique(snps_to_rm))
}

#Identify SNPs that fail HWE
hwe_fails <- filter_hwe_by_pop(geno_data_filt, pop_col = "Superpopulation")
length(hwe_fails) #Remove 679 SNPs

#Remove the SNPs
geno_data_filt <- geno_data_filt %>% 
  select(-all_of(hwe_fails))

# Check how many SNPs are left
print(paste("Original SNPs:", ncol(genotype_columns)))
print(paste("Informative SNPs remaining:", ncol(geno_data_filt %>% select(-SampleID, -Sex, -Population, -Superpopulation))))

#Filter for geographic locations that will likely show differences in lactose persistence based on literature review
geno_matrix_small <- geno_data_filt %>% 
  filter(Superpopulation %in% c("AFR", "EUR", "EAS")) %>% 
  select(-SampleID, -Sex, -Population, -Superpopulation)

#Labels for metadata
labels_small <- geno_data_filt %>% 
  filter(Superpopulation %in% c("AFR", "EUR", "EAS")) %>% 
  select(SampleID, Sex, Population, Superpopulation)

#Calculate the sum of alleles for each SNP
# A sum of 0 means the variant (minor allele) is not present in the 3,202 samples
snp_sums2 <- colSums(geno_matrix_small, na.rm = TRUE)

#Identify SNPs with at least one variant (sum > 0)
#Remove SNPs where everyone has the variant (sum == 2 * number of samples (3202))
#Remove SNPs that are rare in the population (<0.001)
keep_snps2 <- snp_sums2 > 0 & snp_sums2 < 2*nrow(geno_matrix_small) & snp_sums2 > 0.001*nrow(geno_matrix_small)

#Subset data to keep only informative SNPs
genotypes_small_filt <- geno_matrix_small[, names(which(keep_snps2))]
print(paste("Informative SNPs remaining:", ncol(genotypes_small_filt))) #2086

#Check Missing values
sum(is.na(genotypes_small_filt)) #0

rownames(genotypes_small_filt) <- labels_small$SampleID
#PCA =====
#Prepare numeric genotypes (Isolate only the SNP columns for calcs)
geno_matrix <- genotypes_small_filt %>% 
  as.matrix()

#Normalization
#Using scale. = TRUE in prcomp to standardize data
#Scales each SNP to have unit variance (mean = 0, SD = 1)

#Perform PCA
pca_result <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)
#Summary results of PCA
summary_pca <- summary(pca_result)

#Variance explanation (Scree plot)
pca_var <- pca_result$sdev^2
var_explained <- pca_var / sum(pca_var)
cumu_var <- cumsum(var_explained)

#Figures for PCA: Scree plot
plot(var_explained[1:10], type = "b",
     main = "Figure 1: Scree Plot of LCT and MCM6 Region PCA",
     xlab = "Principal Component",
     ylab = "Proportion of Variance",
     col = "blue", pch = 19)

pca_df <- as.data.frame(pca_result$x)
pca_df$Superpopulation <- labels_small$Superpopulation

#Identifying top SNPs ====
#Goal is to identify lactose persistence across genetic groups so want to colour PCA by superpopulation

#PCA of lactose persistence region
ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal()+
  labs(title = "Figure 2: PCA of Lactose Persistence by Region",
       subtitle = "Samples clustered by superpopulation",
       x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"))


#Identify top contributing SNPs (Loadings)
snp_loadings <- data.frame()

#install.packages(c("ggrepel", "factoextra", "caret", "RColorBrewer"))

library(ggrepel)      # Adds non-overlapping text labels to ggplots
library(factoextra)   # Makes nicer dendrograms with fviz_dend()
library(caret)        # Gives us confusionMatrix() and stratified splitting
library(RColorBrewer) # Clean colour palettes

superpop <- labels_small$Superpopulation
print (superpop)

pop <- labels_small$Population
print (pop)

snp_cols <- colnames(geno_matrix)

#give the different super population groups colours
superpop_levels  <- sort(unique(superpop))   # "AFR" "EAS" "EUR"
superpop_colours <- setNames(
  brewer.pal(n = 3, name = "Set1"),
  superpop_levels
)

snp_loadings <- as.data.frame(pca_result$rotation[, 1:3])
snp_loadings$SNP <- rownames(snp_loadings)

snp_loadings$AbsPC1 <- abs(snp_loadings$PC1)
snp_loadings$AbsPC2 <- abs(snp_loadings$PC2)

top_snps_pc1 <- snp_loadings[order(-snp_loadings$AbsPC1), ][1:20, ]

ggplot(top_snps_pc1, aes(x = reorder(SNP, AbsPC1), y = AbsPC1)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  coord_flip() +
  # coord_flip() rotates the chart so SNP names are on the y-axis
  # and are readable — without this they'd overlap on the x-axis
  labs(
    title    = "Figure 3: Top 20 SNPs by Absolute Loading on PC1",
    subtitle = "SNPs with large loadings drive population separation along PC1",
    x        = "SNP",
    y        = "Absolute Loading on PC1"
  ) +
  theme_classic(base_size = 12)

ggsave("fig3_snp_loadings_PC1.png", width = 7, height = 6, dpi = 300)

scree_df <- data.frame(
  PC          = 1:20,
  VarExplPct  = round(var_explained[1:20] * 100, 2),
  CumuVarPct  = round(cumu_var[1:20] * 100, 2)
)

ggplot(scree_df, aes(x = PC, y = VarExplPct)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_line(aes(group = 1), colour = "black", linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:20) +
  labs(
    title    = "Figure 4: Scree Plot — LCT/MCM6 Region PCA",
    subtitle = paste0("PC1 explains ", scree_df$VarExplPct[1],
                      "% of variance; PC2 explains ", scree_df$VarExplPct[2], "%"),
    x        = "Principal Component",
    y        = "% Variance Explained"
  ) +
  theme_classic(base_size = 13)

ggsave("fig4_scree_plot.png", width = 7, height = 4, dpi = 300)

n_pcs_80 <- which(cumu_var >= 0.80)[1]

pca_df$Population <- labels_small$Population


ggplot(pca_df, aes(x = PC1, y = PC3, colour = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_colour_manual(values = superpop_colours) +
  labs(
    title    = "Figure 5: PCA Biplot — PC1 vs PC3",
    subtitle = "Colours = true superpopulation (not used during PCA)",
    x        = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
    y        = paste0("PC3 (", round(var_explained[3]*100, 1), "%)")
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave("fig5_pca_PC1_PC3.png", width = 7, height = 5.5, dpi = 300)


n_pcs_cluster <- min(n_pcs_80, 50)
pca_scores_for_clustering <- pca_result$x[, 1:n_pcs_cluster]

set.seed(42)
sample_idx <- unlist(tapply(
  seq_len(nrow(geno_matrix)),   # indices 1 to however many samples)
  superpop,                      # group each index by its superpopulation
  function(i) sample(i, min(100, length(i)))  # sample up to 100 from each group
))

pca_sub    <- pca_scores_for_clustering[sample_idx, ]
labels_sub <- superpop[sample_idx]
print(table(labels_sub))


dist_mat <- dist(pca_sub, method = "euclidean")
hclust_result <- hclust(dist_mat, method = "ward.D2")

fviz_dend(
  hclust_result,
  k            = 3,           # Cut tree into 3 clusters (one per superpop)
  k_colors     = unname(brewer.pal(3, "Set2")),  # Branch colours
  show_labels  = FALSE,
  rect         = TRUE,        # Draw rectangles around each cluster
  rect_fill    = TRUE,
  main         = "Figure 6: Hierarchical Clustering Dendrogram (Ward's D2)",
  sub          = "k = 3 clusters; n = 300 subsampled individuals"
) +
  theme_classic(base_size = 12)

ggsave("fig6_dendrogram.png", width = 10, height = 6, dpi = 300)


cluster_labels <- cutree(hclust_result, k = 3)

# Cross-tabulation: rows = HC cluster (1/2/3), cols = true superpopulation
cluster_vs_pop <- table(
  Cluster         = cluster_labels,
  Superpopulation = labels_sub
)

print(cluster_vs_pop)
pca_sub_df <- as.data.frame(pca_sub)
pca_sub_df$Superpopulation <- labels_sub
pca_sub_df$Cluster         <- as.factor(cluster_labels)

ggplot(pca_sub_df, aes(x = PC1, y = PC2,
                       colour = Superpopulation,
                       shape  = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = superpop_colours) +
  scale_shape_manual(values = c(16, 17, 15)) +
  # 16 = circle, 17 = triangle, 15 = square — all solid, easy to distinguish
  labs(
    title    = "Figure 7: PCA Scores with HC Cluster Assignments",
    subtitle = "Colour = true superpopulation | Shape = HC cluster (k = 3, Ward's D2)",
    x        = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
    y        = paste0("PC2 (", round(var_explained[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave("fig7_pca_with_clusters.png", width = 8, height = 5.5, dpi = 300)

pca_summary <- pca_df %>%
  group_by(Superpopulation) %>%
  summarise(
    N        = n(),
    Mean_PC1 = round(mean(PC1), 3),
    Mean_PC2 = round(mean(PC2), 3),
    SD_PC1   = round(sd(PC1), 3),
    .groups  = "drop"
  )
print(pca_summary)



#RF 


# Random Forest ====
# SNP column names like "2:135700020:A:G" contain colons and start with a number,
# which breaks R's formula interface. Rename to SNP_1, SNP_2, etc. and store
# the originals so we can restore them on the importance plots
original_snp_names <- colnames(geno_matrix)
safe_snp_names     <- paste0("SNP_", seq_along(original_snp_names))
colnames(geno_matrix) <- safe_snp_names
snp_name_lookup <- setNames(original_snp_names, safe_snp_names)
snp_cols <- safe_snp_names

rf_df <- as.data.frame(geno_matrix)
rf_df$Superpopulation <- as.factor(superpop)

set.seed(123)
train_idx <- createDataPartition(
  rf_df$Superpopulation,
  p    = 0.8,
  list = FALSE
)

train_df <- rf_df[train_idx,  ]
test_df  <- rf_df[-train_idx, ]

print(table(train_df$Superpopulation))
print(table(test_df$Superpopulation))

set.seed(42)
rf_model <- randomForest(
  Superpopulation ~ .,   
  data       = train_df,
  ntree      = 500,
  importance = TRUE,
  do.trace   = 100
)

print(rf_model)

oob_df <- data.frame(
  Trees = seq_len(nrow(rf_model$err.rate)),
  OOB   = rf_model$err.rate[, "OOB"],
  AFR   = rf_model$err.rate[, "AFR"],
  EAS   = rf_model$err.rate[, "EAS"],
  EUR   = rf_model$err.rate[, "EUR"]
)
oob_long <- tidyr::pivot_longer(
  oob_df,
  cols      = -Trees,
  names_to  = "Class",
  values_to = "ErrorRate"
)

ggplot(oob_long, aes(x = Trees, y = ErrorRate, colour = Class)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c("OOB" = "black", superpop_colours)
  ) +
  labs(
    title    = "Figure 8: Random Forest OOB Error vs. Number of Trees",
    subtitle = "Overall OOB error (black) + per-class error rates",
    x        = "Number of Trees",
    y        = "OOB Error Rate",
    colour   = "Class"
  ) +
  theme_classic(base_size = 13)

ggsave("fig8_oob_error.png", width = 7, height = 4.5, dpi = 300)

#tuning with mtry 
cat("\nTuning mtry...\n")
set.seed(42)
tuned_mtry <- tuneRF(
  x          = train_df[, snp_cols],
  y          = train_df$Superpopulation,
  ntreeTry   = 300,
  stepFactor = 2,
  improve    = 0.01,
  trace      = TRUE,
  plot       = TRUE
)

best_mtry <- tuned_mtry[which.min(tuned_mtry[, 2]), 1]
cat("Best mtry:", best_mtry, "\n")

#tunes RF 
cat("\nTraining final Random Forest (tuned mtry =", best_mtry, ")...\n")
set.seed(42)
rf_final <- randomForest(
  Superpopulation ~ .,
  data       = train_df,
  ntree      = 500,
  mtry       = best_mtry,
  importance = TRUE
)
print(rf_final)

#model eval
test_preds <- predict(rf_final, newdata = test_df[, snp_cols])

conf_mat <- confusionMatrix(test_preds, test_df$Superpopulation)
print(conf_mat)

overall_acc   <- round(conf_mat$overall["Accuracy"] * 100, 2)
overall_kappa <- round(conf_mat$overall["Kappa"], 4)

cat("\n=== FINAL TEST SET PERFORMANCE ===\n")
cat("Overall Accuracy:", overall_acc, "%\n")
cat("Cohen's Kappa:   ", overall_kappa, "\n")
cat("(Kappa > 0.8 = excellent; > 0.6 = good; < 0.4 = poor)\n\n")

per_class <- as.data.frame(conf_mat$byClass)[,
                                             c("Sensitivity", "Specificity", "Precision", "F1")]
print(round(per_class, 4))

#conf matrix heatmap 
conf_table        <- as.data.frame(conf_mat$table)
names(conf_table) <- c("Predicted", "Reference", "Freq")

ggplot(conf_table, aes(x = Reference, y = Predicted, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = Freq), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title    = "Figure 9: Confusion Matrix — Random Forest (Test Set)",
    subtitle = paste0("Overall Accuracy: ", overall_acc,
                      "%  |  Cohen's Kappa: ", overall_kappa),
    x        = "True Superpopulation",
    y        = "Predicted Superpopulation",
    fill     = "Count"
  ) +
  theme_classic(base_size = 14)

ggsave("fig9_confusion_matrix.png", width = 6, height = 5, dpi = 300)

#most important SNPs (MDA and MDG (gini)
importance_df     <- as.data.frame(importance(rf_final))
importance_df$SNP <- rownames(importance_df)

# Restore original SNP names for axis labels
importance_df$SNP_label <- snp_name_lookup[importance_df$SNP]

#top 20 by MDA (mean decrease accuracy)
top20_mda <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

ggplot(top20_mda, aes(x = reorder(SNP_label, MeanDecreaseAccuracy),
                      y = MeanDecreaseAccuracy)) +
  geom_col(fill = "darkorange", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Figure 10: Top 20 SNPs — Mean Decrease in Accuracy",
    subtitle = "Higher value = more important for classification",
    x        = "SNP",
    y        = "Mean Decrease in Accuracy"
  ) +
  theme_classic(base_size = 12)

ggsave("fig10_importance_MDA.png", width = 7, height = 6, dpi = 300)

#top 20 by MDG (mean decrease gini)
top20_mdg <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(20)

ggplot(top20_mdg, aes(x = reorder(SNP_label, MeanDecreaseGini),
                      y = MeanDecreaseGini)) +
  geom_col(fill = "mediumpurple", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Figure 11: Top 20 SNPs — Mean Decrease in Gini Impurity",
    subtitle = "Higher value = more splits in all trees relied on this SNP",
    x        = "SNP",
    y        = "Mean Decrease in Gini"
  ) +
  theme_classic(base_size = 12)

ggsave("fig11_importance_Gini.png", width = 7, height = 6, dpi = 300)

#top 10 SNPs by MDA
print(top20_mda[1:10, c("SNP_label", "MeanDecreaseAccuracy", "MeanDecreaseGini")])
