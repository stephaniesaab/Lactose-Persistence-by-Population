#Stats assignment 3 sh

#Chromosome region:
#2:135,700,000-136,000,000
# This downloads ONLY the requested 300,000 base pairs from the IGSR server
tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz 2:135700000-136000000 > lactose_region.vcf

#Metadata (population, sex)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt


# 1. Download the VCF index first (required for the slice)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

# 2. Use tabix to download the Lactose Persistence chromosome region only (approx 300kb)
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 2:135700000-136000000 > lactose_data.vcf
