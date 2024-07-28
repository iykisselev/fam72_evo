#!/usr/bin/env bash
# Trim and split chimpanzee reads

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
metadata="$repo_root/data/metadata/chimpanzee_samples.tsv"
output_dir="$repo_root/data/derived/3.1/chimpanzee"

mkdir -p "$output_dir/trimmed" "$output_dir/chunks"

tail -n +2 "$metadata" | while IFS=$'\t' read -r sample_id subspecies fastq_r1 fastq_r2; do
  trim_r1="$output_dir/trimmed/${sample_id}_R1.trimmed.fastq.gz"
  trim_r2="$output_dir/trimmed/${sample_id}_R2.trimmed.fastq.gz"

  trimmomatic PE -phred33 -threads 8 \
    "$fastq_r1" "$fastq_r2" \
    "$trim_r1" "$output_dir/trimmed/${sample_id}_R1.unpaired.fastq.gz" \
    "$trim_r2" "$output_dir/trimmed/${sample_id}_R2.unpaired.fastq.gz" \
    SLIDINGWINDOW:4:30 MINLEN:36

  sample_chunk_dir="$output_dir/chunks/$sample_id"
  mkdir -p "$sample_chunk_dir/r1" "$sample_chunk_dir/r2"
  seqkit split2 --by-part 20 --out-dir "$sample_chunk_dir/r1" "$trim_r1"
  seqkit split2 --by-part 20 --out-dir "$sample_chunk_dir/r2" "$trim_r2"
done
