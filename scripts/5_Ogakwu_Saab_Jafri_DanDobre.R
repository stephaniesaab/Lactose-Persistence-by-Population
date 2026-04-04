#Assignment 3 for Statistics for Bioinformatics: Gene Classifications

#Stats assignment 3 sh
#This is the bash commands used to download the data
#Chromosome region:
#2:135,700,000-136,000,000
# This downloads only the requested 300,000 base pairs from the IGSR server
#tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr2:135700000-136000000 > lactose_region.vcf

#Metadata (population, sex)
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

####################################################################################
#Download
#wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

#Rename
#mv 1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz chr2_full.vcf.gz

#Create .tbi index
#tabix -p vcf chr2_full.vcf.gz

#Get small slice for lactose gene
#tabix -h chr2_full.vcf.gz chr2:135700000-136000000 > lactose_persistence.vcf

#Or do it from URL:
#Slice the lactose region directly from the remote URL (This gave a file with 111 lines of comments)
#tabix -h https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV_phased_panel.vcf.gz chr2:135700000-136000000 > lactose_persistence_slice.vcf

#Check if there is data:
#grep -v "^#" lactose_persistence.vcf | wc -l #Success: got 6776 lines


#Load libraries ====
library(vcfR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(randomForest)
library(pegas)
library(caret)
library(RColorBrewer)
library(plotly)
library(factoextra)

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
#Should output processed variant: ~26319

data_dir <- "../data"
figures_dir <- "../figures"
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

vcf_candidates <- c(
  file.path(data_dir, "lactose_persistence_slice_large.vcf"),
  file.path(data_dir, "lactose_persistence.vcf"),
  file.path(data_dir, "lactose_persistence_slice.vcf")
)

vcf_path <- vcf_candidates[file.exists(vcf_candidates)][1]

if (length(vcf_path) == 0 || is.na(vcf_path)) {
  stop("No VCF file found. Expected one of: lactose_persistence_slice_large.vcf, lactose_persistence.vcf, lactose_persistence_slice.vcf in ../data/")
}

cat("Using VCF file:", vcf_path, "\n")
vcf <- read.vcfR(vcf_path, verbose = TRUE)

#Read metadata
metadata <- read.table(file.path(data_dir, "20130606_g1k_3202_samples_ped_population.txt"), header = TRUE, sep = " ") #Columns separated by spaces
names(metadata)

#Extract genotypes and sample IDs
#Extract genotype matrix (Converts genotypes into number of alternate alleles)
genotypes <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
nrow(genotypes)

#Transpose so samples are rows and SNPs are columns
genotypes_t <- t(genotypes)
ncol(genotypes_t)

#Convert data to a dataframe to make sampleIDs a column
geno_df <- as.data.frame(genotypes_t) #nrow  = 
geno_df$SampleID <- rownames(geno_df)
geno_df <- geno_df %>%
  filter(SampleID %in% metadata$SampleID)

#Merge metadata with genotypes
geno_final <- geno_df %>% 
  inner_join(metadata, by = c("SampleID" = "SampleID"))
nrow(geno_final)
ncol(geno_final)

#Check populations
print(table(geno_final$Population)) #Three letter codes for locations
print(table(geno_final$Superpopulation)) #Africa, America, East Asia, Europe, South Asia


#Write out CSV with variants and data
#write.csv(geno_final, "../data/genotype_data.csv")

# EDA ====
#Exploratory data analysis to identify which subpopulations to extract before doing clustering by genotypes
##1. QC ====
# Population Sample Size Bar Chart
pop_counts <- geno_final %>%
  count(Superpopulation, Population)

p_sample_pop <- ggplot(pop_counts, aes(x = reorder(Population, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Sample Size by Population",
    x = "Population",
    y = "Sample Count",
    fill = "Superpopulation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_sample_pop)
ggsave(file.path(figures_dir, "fig_sample_size_pop.png"), plot = p_sample_pop, width = 9, height = 5, dpi = 300)

# Population Counts by Superpopulation
p_sample_superpop <- ggplot(pop_counts %>% group_by(Superpopulation) %>% summarise(n = sum(n)),
       aes(x = reorder(Superpopulation, -n), y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  labs(title = "Sample Size by Superpopulation", x = "Superpopulation", y = "Count") +
  theme_minimal()
print(p_sample_superpop)
ggsave(file.path(figures_dir, "fig_sample_size_superpop.png"), plot = p_sample_superpop, width = 6, height = 4.5, dpi = 300)

## Per-Sample Heterozygosity
het_rates <- rowMeans(geno_final == 1)
het_df <- data.frame(
  SampleID        = geno_final$SampleID,
  Heterozygosity  = het_rates,
  Superpopulation = geno_final$Superpopulation
)

het_mean <- mean(het_df$Heterozygosity)
het_sd   <- sd(het_df$Heterozygosity)
het_df$Outlier <- abs(het_df$Heterozygosity - het_mean) > 3 * het_sd

cat("Heterozygosity outliers:", sum(het_df$Outlier), "\n")
print(het_df[het_df$Outlier, ])

p_het <- ggplot(het_df, aes(x = Superpopulation, y = Heterozygosity, fill = Superpopulation)) +
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
print(p_het)
ggsave(file.path(figures_dir, "fig_heterozygosity.png"), plot = p_het, width = 8, height = 5, dpi = 300)


##2. Variant Density Across the Region ====
snp_positions <- data.frame(
  CHROM = getCHROM(vcf),
  POS   = as.numeric(getPOS(vcf)),
  ID    = getID(vcf)
)

#Define slice coordinates
region_start <- min(snp_positions$POS, na.rm = TRUE)
region_end   <- max(snp_positions$POS, na.rm = TRUE)

#Filter the SNP positions to only include this slice
snp_positions_slice <- snp_positions %>%
  filter(POS >= region_start & POS <= region_end)

# Plot the specific region
p_variant_density <- ggplot(snp_positions_slice, aes(x = POS)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  labs(
    title = "Figure 3: Variant Density Across LCT/MCM6 Region",
    subtitle = paste0("Genomic Slice: ", scales::comma(region_start), " - ", scales::comma(region_end), " bp"),
    x = "Genomic Position (Chromosome 2)",
    y = "Number of Variants"
  ) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal()
print(p_variant_density)
ggsave(file.path(figures_dir, "fig3_variant_density.png"), plot = p_variant_density, width = 8, height = 4.5, dpi = 300)

#Initial Subsetting and cleaning ====

#Filter for target populations first (AFR, EUR, EAS, SAS).
#Ensures HWE and MAF calcs are for specific analysis groups.

geno_subset <- geno_final %>% 
  filter(Superpopulation %in% c("AFR", "EUR", "EAS"))

#Separate genotype columns from the metadata columns
metadata_cols <- c("SampleID", "FamilyID", "MotherID", "FatherID", "Sex", "Population", "Superpopulation")
all_snps <- setdiff(names(geno_subset), metadata_cols)
length(all_snps)

##1. SNP Filtering (MAF and informative) ====
#Remove rare variants and monomorphic loci
#Calculate the sum of alleles for each SNP
# A sum of 0 means the variant (minor allele) is not present in the 3,202 samples
mat_temp <- as.matrix(geno_subset[, all_snps])
snp_sums <- colSums(mat_temp)
n_samples <- nrow(geno_subset)
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
length(active_snps)
geno_matrix <- as.matrix(geno_final[, active_snps])
rownames(geno_matrix) <- geno_final$SampleID
head(geno_matrix[, 1:5])

print(paste("Final SNPs for analysis:", length(active_snps)))

##2. HWE filtering ====
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
# PCA =====
#Prepare numeric genotypes (Isolate only the SNP columns for calcs)

#Normalization
#Using scale. = TRUE in prcomp to standardize data
#Scales each SNP to have unit variance (mean = 0, SD = 1)

#Perform PCA
pca_result <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)

#Summary results of PCA
summary_pca <- summary(pca_result)
print(summary_pca)

#Variance explanation (Scree plot)
pca_var <- pca_result$sdev^2
var_explained <- pca_var / sum(pca_var)
cumu_var <- cumsum(var_explained)

cat("Variance explained by PC1:", round(var_explained[1] * 100, 2), "%\n")
cat("Variance explained by PC2:", round(var_explained[2] * 100, 2), "%\n")
cat("Variance explained by PC3:", round(var_explained[3] * 100, 2), "%\n")
cat("Cumulative variance (PC1-5):", round(cumu_var[5] * 100, 2), "%\n")

#Figures for PCA: Scree plot
png(file.path(figures_dir, "fig1_scree_base.png"), width = 900, height = 600, res = 120)
plot(var_explained[1:10], type = "b",
     main = "Figure 1: Scree Plot of LCT and MCM6 Region PCA",
     xlab = "Principal Component",
     ylab = "Proportion of Variance",
     col = "blue", pch = 19)
dev.off()

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

p_pca_pre_pc1pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA of Genotype Data (pre-HWE filtering)",
    subtitle = "Colored by Superpopulation",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC2 (", pca_var_pct[2], "% variance)")
  ) +
  theme_minimal()
print(p_pca_pre_pc1pc2)
ggsave(file.path(figures_dir, "fig_pca_preHWE_PC1PC2.png"), plot = p_pca_pre_pc1pc2, width = 7, height = 5, dpi = 300)

p_pca_pre_pc1pc3 <- ggplot(pca_df, aes(x = PC1, y = PC3, color = Superpopulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "PCA: PC1 vs PC3 (pre-HWE filtering)",
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC3 (", pca_var_pct[3], "% variance)")
  ) +
  theme_minimal()
print(p_pca_pre_pc1pc3)
ggsave(file.path(figures_dir, "fig_pca_preHWE_PC1PC3.png"), plot = p_pca_pre_pc1pc3, width = 7, height = 5, dpi = 300)

scree_df <- data.frame(PC = 1:min(20, length(pca_var_pct)), VariancePct = pca_var_pct[1:min(20, length(pca_var_pct))])
p_scree <- ggplot(scree_df, aes(x = PC, y = VariancePct)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot (pre-HWE filtering)", x = "Principal Component", y = "Variance Explained (%)") +
  scale_x_continuous(breaks = scree_df$PC) +
  theme_minimal()
print(p_scree)
ggsave(file.path(figures_dir, "fig_scree_ggplot.png"), plot = p_scree, width = 8, height = 4.5, dpi = 300)

## 1.Identifying top SNPs ====
#Goal is to identify lactose persistence across genetic groups so want to colour PCA by superpopulation

#PCA of lactose persistence region
p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal()+
  labs(title = "Figure 2: PCA of Lactose Persistence by Region",
       subtitle = "Samples clustered by superpopulation",
       x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"))
print(p_pca_main)
ggsave(file.path(figures_dir, "fig2_pca_main.png"), plot = p_pca_main, width = 7, height = 5, dpi = 300)


## 1.Identify top contributing SNPs (Loadings) ====
snp_loadings <- data.frame(
  SNP = rownames(pca_result$rotation),
  PC1_loading = abs(pca_result$rotation[,1])
)

top_snps <- snp_loadings %>% 
  arrange(desc(PC1_loading)) %>% 
  head(10)
#Get top 10 SNPs
print(top_snps)



#K-means ====
set.seed(1516)

# Take only the first 10 PCs to capture the main signal
# Take only the first 11 PCs to capture the main signal
n_pcs_for_km <- min(11, ncol(pca_result$x))
pca_for_km <- pca_result$x[, 1:n_pcs_for_km, drop = FALSE]

#Get optimal clusters

# Calculate WSS for K=1 to K=10
#Expected 3 or 5 based on Superpopulations (AFR, EUR, EAS)
wss <- sapply(1:10, function(k){
  kmeans(pca_for_km, centers = k, nstart = 20, iter.max = 100)$tot.withinss
})

# The Elbow Plot
elbow_df <- data.frame(K = 1:10, WSS = wss)
p_elbow <- ggplot(elbow_df, aes(x = K, y = WSS)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Figure 5: Elbow Method for Optimal K",
       x = "Number of Clusters (K)", y = "Total Within-Cluster SS") +
  theme_minimal()
print(p_elbow)
ggsave(file.path(figures_dir, "fig5_elbow.png"), plot = p_elbow, width = 7, height = 5, dpi = 300)

#Elbow plot had no clear Bend, use silhouette method to get number of clusters to use
p_silhouette <- fviz_nbclust(pca_for_km, kmeans, method = "silhouette") +
  labs(title = "Optimal Number of Clusters (Silhouette Method)")
print(p_silhouette)
ggsave(file.path(figures_dir, "fig6_silhouette.png"), plot = p_silhouette, width = 7, height = 5, dpi = 300)

#Indicated around 9 clusters, but showed average silhouette increased greatly around 6

#Run kmeans on SNP columns, 6 clusters based on elbow plot
#Run kmeans on SNP columns, 6 clusters based on elbow plot and Silhouette
km_res <- kmeans(pca_for_km, centers = 6, nstart = 25)

#Add cluster assignments for plotting
geno_final$Cluster <- as.factor(km_res$cluster)

#Check results against the actual populations
print(table(geno_final$Superpopulation, geno_final$Cluster))

#RF ====

# SNP column names like "2:135700020:A:G" contain colons and start with a number,
# which breaks R's formula interface. Rename to SNP_1, SNP_2, etc. and store
# the originals so we can restore them on the importance plots
superpop <- geno_final$Superpopulation
superpop_levels  <- sort(unique(superpop))
superpop_colours <- setNames(
  brewer.pal(n = length(superpop_levels), name = "Set1"),
  superpop_levels
)

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

p_oob <- ggplot(oob_long, aes(x = Trees, y = ErrorRate, colour = Class)) +
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
print(p_oob)
ggsave(file.path(figures_dir, "fig8_oob_error.png"), plot = p_oob, width = 7, height = 4.5, dpi = 300)

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

p_conf <- ggplot(conf_table, aes(x = Reference, y = Predicted, fill = Freq)) +
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
print(p_conf)
ggsave(file.path(figures_dir, "fig9_confusion_matrix.png"), plot = p_conf, width = 6, height = 5, dpi = 300)

#most important SNPs (MDA and MDG (gini)
importance_df     <- as.data.frame(importance(rf_final))
importance_df$SNP <- rownames(importance_df)

# Restore original SNP names for axis labels
importance_df$SNP_label <- snp_name_lookup[importance_df$SNP]

#top 20 by MDA (mean decrease accuracy)
top20_mda <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

p_mda <- ggplot(top20_mda, aes(x = reorder(SNP_label, MeanDecreaseAccuracy),
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
print(p_mda)
ggsave(file.path(figures_dir, "fig10_importance_MDA.png"), plot = p_mda, width = 7, height = 6, dpi = 300)

#top 20 by MDG (mean decrease gini)
top20_mdg <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(20)

p_mdg <- ggplot(top20_mdg, aes(x = reorder(SNP_label, MeanDecreaseGini),
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
print(p_mdg)
ggsave(file.path(figures_dir, "fig11_importance_Gini.png"), plot = p_mdg, width = 7, height = 6, dpi = 300)

#top 10 SNPs by MDA
print(top20_mda[1:10, c("SNP_label", "MeanDecreaseAccuracy", "MeanDecreaseGini")])

cat("\nDone. All figures saved to ../figures/\n")
