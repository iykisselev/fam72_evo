#!/usr/bin/env bash
# Chapter 3.1
# Generate the BEAST alignment and run configured BEAST analyses

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
reference="$repo_root/data/raw/reference/GRCh38.primary_assembly.genome.fa"
gene_manifest="$repo_root/data/metadata/fam72.grch38.tsv"
vcf_dir="$repo_root/data/derived/3.2/loci"
beast_dir="$repo_root/data/derived/3.1/beast"
chimp_fasta="$repo_root/data/derived/3.1/chimpanzee/consensus/predominant_fam72a_haplotypes.fa"
gorilla_fasta="$repo_root/data/raw/gorilla/Gorilla_gorilla_FAM72A.fa"
result_dir="$repo_root/results/3.1/beast"

mkdir -p "$beast_dir/consensus" "$result_dir"
human_fasta="$beast_dir/consensus/human_haplotypes_above_100.fa"
: > "$human_fasta"

tail -n +2 "$beast_dir/human_haplotypes_above_100.tsv" | while IFS=$'\t' read -r gene haplotype count sample copy; do
  vcf="$vcf_dir/${gene}.1kg.vcf.gz"
  region="$(awk -F $'\t' -v target="$gene" '$4 == target {print $1 ":" $2 "-" $3}' "$gene_manifest")"
  printf '>%s.%s\n' "$(tr '[:upper:]' '[:lower:]' <<< "$gene")" "$haplotype"
  bcftools consensus --fasta-ref "$reference" --regions "$region" --sample "$sample" --haplotype "$copy" "$vcf" | sed '1d'
  >> "$human_fasta"
done

combined_fasta="$beast_dir/fam72_human_chimp_gorilla_top_haplotypes.fa"
cat "$human_fasta" "$chimp_fasta" "$gorilla_fasta" > "$combined_fasta"
mafft --localpair --maxiterate 1000 "$combined_fasta" > "$beast_dir/fam72_human_chimp_gorilla_top_haplotypes.aligned.fa"
iqtree2 --sequence "$beast_dir/fam72_human_chimp_gorilla_top_haplotypes.aligned.fa" --model MFP --redo

# Run nested-sampling XMLs first. The XMLs must specify 4 particles, chain length
# 10^5, subchain length 2x10^4, stated calibrations, and one strict/one relaxed clock.
for xml in "$repo_root"/data/metadata/beast/*nested_sampling.xml; do
  [[ -f "$xml" ]] || continue
  beast -working "$result_dir" "$xml"
done

# Run the selected clock model three times after examining nested-sampling results.
for xml in "$repo_root"/data/metadata/beast/*selected_clock_run*.xml; do
  [[ -f "$xml" ]] || continue
  beast -working "$result_dir" "$xml"
done
