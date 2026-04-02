#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Load libraries ====
library(vcfR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(randomForest)
library(pegas)

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
