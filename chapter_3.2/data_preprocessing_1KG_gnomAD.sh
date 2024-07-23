#!/bin/bash
# The script extracts and filters FAM72A-D genomic regions from the 1000 Genomes and V3.1.2 gnomAD VCF files based on coordinates provided in 'fam72.grch38.tsv'. It downloads data from the 1000 Genomes Project, then extracts and processes genomic regions to keep only biallelic SNPs and remove parents from trio sample groups. Filtering is not applied to the GNOMAD data.

# Download necessary files from 1000 Genomes Project and gnomAD
download_files() {
  local urls=(
	"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
	"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt"
	"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz.tbi"
  )

  for url in "${urls[@]}"; do
	wget -q "$url"
  done
}

# Extract non-zero parent IDs to remove them from future analyses
extract_parents_to_remove() {
  awk 'NR > 1 {if ($2 != 0 && $2 != "") print $2; if ($3 != 0 && $3 != "") print $3;}' 1kGP.3202_samples.pedigree_info.txt | \
  sort | uniq > samples_to_remove.txt
}

# Function to extract and filter VCF for a specific region from 1000 Genomes
extract_1000G_region() {
  local region="$1:$2-$3"
  local output_file="$4.1KG.vcf.gz"

  tabix -h "$vcf_1KG" "$region" | \
  bcftools norm -m +any -Ou | \
  bcftools view -S ^samples_to_remove.txt --force-samples -m2 -M2 -v snps -Oz -o "$output_file"
}

# Function to extract VCF for a specific region from gnomAD
extract_gnomad_region() {
  local region="$1:$2-$3"
  local output_file="$4.gnomad.vcf.gz"

  tabix -h "$vcf_gnomad" "$region" | bgzip -c > "$output_file"
}

# Main function to process regions
process_regions() {
  while read -r chrom start end name; do
	extract_1000G_region "$chrom" "$start" "$end" "$name"
	extract_gnomad_region "$chrom" "$start" "$end" "$name"
  done < fam72.grch38.tsv
}

# Load genomic tools
module load tabix bcftools

# Define path to the VCF files
vcf_1KG="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"
vcf_gnomad="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz"

# Execute functions
download_files
extract_parents_to_remove
process_regions
