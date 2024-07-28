#!/usr/bin/env bash
# Chapter 3.1
# Varinat calling with FreeBayes. Joint variant call on chromosome 1 separately
# for P. t. troglodytes, P. t. schweinfurthii, and P. t. verus.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
metadata="$repo_root/data/metadata/chimpanzee_samples.tsv"
reference="$repo_root/data/raw/chimpanzee/reference/NHGRI_mPanTro3-v1.1.fa"
output_dir="$repo_root/data/derived/3.1/chimpanzee"

mkdir -p "$output_dir/vcf"

for subspecies in troglodytes schweinfurthii verus; do
  mapfile -t bams < <(
    awk -F $'\t' -v group="$subspecies" 'NR > 1 && $2 == group {print $1}' "$metadata" \
      | sed "s#^#$output_dir/bam/#; s#\$#.chr1.filtered.bam#"
  )

  bam_list="$output_dir/bam/${subspecies}.bam.list"
  raw_vcf="$output_dir/vcf/${subspecies}.raw.vcf.gz"
  filtered_vcf="$output_dir/vcf/${subspecies}.filtered.vcf.gz"
  printf '%s\n' "${bams[@]}" > "$bam_list"

  freebayes --fasta-reference "$reference" --region chr1 --bam-list "$bam_list" \
    | bgzip --stdout > "$raw_vcf"
  tabix --force --preset vcf "$raw_vcf"

  bcftools view \
    --include 'QUAL > 19 && QUAL / INFO/AO > 10 && INFO/SAF > 0 && INFO/SAR > 0 && INFO/RPR > 1 && INFO/RPL > 1' \
    --output-type u "$raw_vcf" \
    | bcftools norm --multiallelics -any --output-type u \
    | bcftools view --genotype ^miss --output-type z --output "$filtered_vcf"
  tabix --force --preset vcf "$filtered_vcf"
done
