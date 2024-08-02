#!/usr/bin/env bash
# Chapter 3.3
# For each population, an artificial homozygous ancestral sample is appended from INFO/AA. 
# glactools then converts the VCF to ACF, marks that sample as both ancestor and root, and exports the unfolded BetaScan format.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_dir="$repo_root/data/derived/3.3/ihs/populations"
reference_fai="$repo_root/data/raw/reference/GRCh38.primary_assembly.genome.fa.fai"
output_dir="$repo_root/data/derived/3.3/betascan"


mkdir -p "$output_dir"/{dummy_vcfs,acf,betascan_input}

for input_vcf in "$input_dir"/*.ihs_input.vcf.gz; do
  [[ -s "$input_vcf" ]] || continue
  population="$(basename "$input_vcf" .ihs_input.vcf.gz)"
  dummy_vcf="$output_dir/dummy_vcfs/${population}.with_ancestral.vcf.gz"

  # Add an ANCESTRAL sample. FORMAT fields other than GT are emitted as missing.
  bcftools view --output-type v "$input_vcf" | awk 'BEGIN {FS=OFS="\t"}
    /^##/ {print; next}
    /^#CHROM/ {print $0, "ANCESTRAL"; next}
    {
      split($8, info_items, ";"); aa=""
      for (i in info_items) if (info_items[i] ~ /^AA=/) {split(info_items[i], item, "="); aa=item[2]}
      sub(/\|.*/, "", aa)
      if (aa == "" || aa == ".") next
      allele_index=-1
      if (aa == $4) allele_index=0
      else {split($5, alternate, ","); for (i in alternate) if (aa == alternate[i]) allele_index=i}
      if (allele_index < 0) next
      split($9, format_fields, ":"); sample=""
      for (i=1; i<=length(format_fields); i++) {
        value=(format_fields[i] == "GT" ? allele_index "/" allele_index : ".")
        sample=sample (i == 1 ? "" : ":") value
      }
      print $0, sample
    }' | bgzip --stdout > "$dummy_vcf"
  tabix --force --preset vcf "$dummy_vcf"

  acf="$output_dir/acf/${population}.acf.gz"
  rooted_acf="$output_dir/acf/${population}.rooted.acf.gz"
  glactools vcfm2acf --onlyGT --fai "$reference_fai" "$dummy_vcf" > "$acf"
  glactools usepopsrootanc --root ANCESTRAL --anc ANCESTRAL "$acf" > "$rooted_acf"
  glactools acf2betascan --useanc "$rooted_acf" | gzip --stdout \
    > "$output_dir/betascan_input/${population}.betascan.txt.gz"
done
