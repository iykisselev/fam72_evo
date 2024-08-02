#!/usr/bin/env bash
# Chapter 3.3
# This script uses released FILTER and INFO fields of the gNOMAD vcf file release

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_vcf="$repo_root/data/raw/gnomad/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz"
manifest="$repo_root/data/metadata/fam72.grch38.tsv"
output_dir="$repo_root/data/derived/3.3/gnomad"

mkdir -p "$output_dir"

while IFS=$'\t' read -r chrom start end gene; do
  [[ "$chrom" == "chrom" ]] && continue
  output_vcf="$output_dir/${gene}.gnomad_v3.1.2.filtered.vcf.gz"

  bcftools view \
    --regions "${chrom}:${start}-${end}" \
    --apply-filters PASS \
    --include 'INFO/AC>0 && INFO/InbreedingCoeff>=-0.3 && ((TYPE="snp" && INFO/AS_VQSLOD>=-2.7739) || (TYPE="indel" && INFO/AS_VQSLOD>=-1.0606))' \
    --output-type z --output "$output_vcf" "$source_vcf"
  tabix --force --preset vcf "$output_vcf"
done < "$manifest"
