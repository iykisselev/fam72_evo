#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_vcf="$repo_root/data/derived/3.2/chr1.1kg.parents_removed.multiallelic_snps.vcf.gz"
ancestral_vcf="$repo_root/data/raw/ensembl/ensembl_110_homo_sapiens_chr1_ancestral.vcf.gz"
metadata="$repo_root/data/metadata/population_info.tsv"
output_dir="$repo_root/data/derived/3.3/ihs"

mkdir -p "$output_dir/populations" "$output_dir/sample_lists"
annotated_vcf="$output_dir/chr1.1kg.ancestral_allele_snps.vcf.gz"

bcftools annotate --annotations "$ancestral_vcf" --columns INFO/AA --output-type u "$input_vcf" \
  | bcftools view --include 'INFO/AA!="." && INFO/AA!=""' --output-type z --output "$annotated_vcf"
tabix --force --preset vcf "$annotated_vcf"

tail -n +2 "$metadata" | cut -f2 | sort -u | while read -r population; do
  sample_list="$output_dir/sample_lists/${population}.txt"
  awk -F $'\t' -v population="$population" 'NR > 1 && $2 == population {print $1}' "$metadata" > "$sample_list"
  bcftools view --samples-file "$sample_list" --output-type z \
    --output "$output_dir/populations/${population}.ihs_input.vcf.gz" "$annotated_vcf"
  tabix --force --preset vcf "$output_dir/populations/${population}.ihs_input.vcf.gz"
done
