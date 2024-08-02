#!/usr/bin/env bash
# Chapter 3.3

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_vcf="$repo_root/data/derived/3.2/chr1.1kg.parents_removed.multiallelic_snps.vcf.gz"
output_dir="$repo_root/data/derived/3.3/neutrality/background_vcfs"
window_manifest="$repo_root/data/derived/3.3/neutrality/background_segments.tsv"
chromosome_length=248956422
segment_width=1000000

mkdir -p "$output_dir"
printf 'chrom\tstart\tend\tsegment\n' > "$window_manifest"

segment=1
for ((start = 1; start <= chromosome_length; start += segment_width)); do
  end=$((start + segment_width - 1))
  (( end > chromosome_length )) && end=$chromosome_length
  name=$(printf 'chr1_%09d_%09d' "$start" "$end")
  printf 'chr1\t%s\t%s\t%s\n' "$start" "$end" "$name" >> "$window_manifest"
  bcftools view --regions "chr1:${start}-${end}" --output-type z \
    --output "$output_dir/${name}.vcf.gz" "$input_vcf"
  tabix --force --preset vcf "$output_dir/${name}.vcf.gz"
  ((segment += 1))
done
