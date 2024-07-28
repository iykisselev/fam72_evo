#!/usr/bin/env bash
# Chapter 3.1
# SHAPEIT phasing of chimp variants calls using subspecies-specific population sizes

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
gene_manifest="$repo_root/data/metadata/fam72.grch38.tsv"
reference="$repo_root/data/raw/chimpanzee/reference/NHGRI_mPanTro3-v1.1.fa"
shapeit_map="$repo_root/data/raw/chimpanzee/recombination/chimpanzee_chr1.map"
output_dir="$repo_root/data/derived/3.1/chimpanzee"

mkdir -p "$output_dir/phased" "$output_dir/consensus"

Rscript "$repo_root/chapter_3.1/02_shift_chimpanzee_multiallelic_positions.R"

for subspecies in troglodytes schweinfurthii; do
  shifted_vcf="$output_dir/vcf/${subspecies}.shifted.vcf.gz"
  if [[ "$subspecies" == "schweinfurthii" ]]; then
    effective_size=32492
  else
    effective_size=44000
  fi

  prefix="$output_dir/phased/${subspecies}"
  shapeit \
    --input-vcf "$shifted_vcf" \
    --input-map "$shapeit_map" \
    --window 0.5 \
    --rho 0.00119 \
    --effective-size "$effective_size" \
    --output-max "$prefix.haps" "$prefix.sample"
  shapeit --convert --input-haps "$prefix.haps" "$prefix.sample" \
    --output-vcf "$prefix.shifted.phased.vcf.gz"
  tabix --force --preset vcf "$prefix.shifted.phased.vcf.gz"
done

Rscript "$repo_root/chapter_3.1/03_restore_chimpanzee_positions.R"

fam72a_region="$(awk -F $'\t' '$4 == "FAM72A" {print $1 ":" ($2 - 5000) "-" ($3 + 5000)}' "$gene_manifest")"
for vcf in "$output_dir"/phased/*.restored.vcf.gz; do
  subspecies="$(basename "$vcf" .restored.vcf.gz)"
  while read -r sample; do
    for haplotype in 1 2; do
      fasta="$output_dir/consensus/${subspecies}.${sample}.h${haplotype}.fa"
      {
        printf '>%s.%s.h%s\n' "$subspecies" "$sample" "$haplotype"
        bcftools consensus --fasta-ref "$reference" --regions "$fam72a_region" \
          --sample "$sample" --haplotype "$haplotype" "$vcf" | sed '1d'
      } > "$fasta"
    done
  done < <(bcftools query --list-samples "$vcf")
done
