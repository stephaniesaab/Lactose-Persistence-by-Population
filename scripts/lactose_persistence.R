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

