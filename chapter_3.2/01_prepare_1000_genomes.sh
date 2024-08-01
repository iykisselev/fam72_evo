#!/usr/bin/env bash
# Chapter 3.2

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="$repo_root/data/raw/1000_genomes"
metadata_dir="$repo_root/data/metadata"
derived_dir="$repo_root/data/derived/3.2"

source_vcf="$raw_dir/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"
pedigree="$raw_dir/1kGP.3202_samples.pedigree_info.txt"
gene_manifest="$metadata_dir/fam72.grch38.tsv"
locus_manifest="$metadata_dir/fam72_srgap2_loci.grch38.tsv"

mkdir -p "$derived_dir/loci"

# The published pedigree table has father and mother identifiers in columns 2 and 3.
# Parents are removed so each trio contributes offspring only
awk 'NR > 1 {if ($2 != "0" && $2 != "") print $2; if ($3 != "0" && $3 != "") print $3}' \
  "$pedigree" | sort -u > "$derived_dir/trio_parents_to_remove.txt"

# Keep SNPs and remove parents before merging biallelic records.
bcftools view \
  --samples-file "^$derived_dir/trio_parents_to_remove.txt" \
  --force-samples \
  --types snps \
  --min-alleles 2 \
  --max-alleles 2 \
  --output-type u \
  "$source_vcf" \
  | bcftools norm --multiallelics +any --output-type z \
      --output "$derived_dir/chr1.1kg.parents_removed.multiallelic_snps.vcf.gz"

tabix --force --preset vcf "$derived_dir/chr1.1kg.parents_removed.multiallelic_snps.vcf.gz"

extract_manifest_regions() {
  local manifest="$1"
  local label_column="$2"

  while IFS=$'\t' read -r chrom start end label; do
    [[ "$chrom" == "chrom" ]] && continue
    [[ -n "$label" ]] || continue

    bcftools view \
      --regions "${chrom}:${start}-${end}" \
      --output-type z \
      --output "$derived_dir/loci/${label}.1kg.vcf.gz" \
      "$derived_dir/chr1.1kg.parents_removed.multiallelic_snps.vcf.gz"
    tabix --force --preset vcf "$derived_dir/loci/${label}.1kg.vcf.gz"
  done < "$manifest"
}

extract_manifest_regions "$gene_manifest" "gene"
extract_manifest_regions "$locus_manifest" "locus"
