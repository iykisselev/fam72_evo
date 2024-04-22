#!/bin/bash
# The code extracts and filters FAM72A-D genomic regions from the 1000 Genomes VCF file based on coordinates provided in 'fam72.grch38.tsv'. It first downloads data from the 1000 Genomes Project, then extracts and processes genomic regions to keep only biallelic SNPs and remove parents from trio sample groups.
# Download index and pedigree and population information files from 1000 Genomes Project
wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

# Extract non-zero parent IDs to remove them from future analyses
awk 'NR > 1 {if ($2 != 0 && $2 != "") print $2; if ($3 != 0 && $3 != "") print $3;}' 1kGP.3202_samples.pedigree_info.txt | \
sort | uniq > samples_to_remove.txt

# Load genomic tools
module load tabix bcftools

# Define path to the VCF file
vcf_1KG="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"

# Function to extract and filter VCF for a specific region
function extract_region() {
    local region="$1:$2-$3"
    local output_file="$4.1KG.vcf.gz"

    # Extract region, normalize multiallelic sites, filter out parents and select only biallelic SNPs
    tabix -h "$vcf_1KG" "$region" | \
    bcftools norm -m +any -Ou | \
    bcftools view -S ^samples_to_remove.txt --force-samples -m2 -M2 -v snps -Oz -o "$output_file"
}

# Loop through each line in the coordinates file
while read -r chrom start end name; do
    extract_region "$chrom" "$start" "$end" "$name"
done < fam72.grch38.tsv
