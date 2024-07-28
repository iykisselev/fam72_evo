#!/usr/bin/env bash
# Map, sort, deduplicate, and filter chimpanzee reads

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
metadata="$repo_root/data/metadata/chimpanzee_samples.tsv"
reference="$repo_root/data/raw/chimpanzee/reference/NHGRI_mPanTro3-v1.1.fa"
output_dir="$repo_root/data/derived/3.1/chimpanzee"

mkdir -p "$output_dir/bam"

tail -n +2 "$metadata" | while IFS=$'\t' read -r sample_id subspecies fastq_r1 fastq_r2; do
  sample_chunk_dir="$output_dir/chunks/$sample_id"
  chunk_bams=()

  for r1_chunk in "$sample_chunk_dir"/r1/*; do
    chunk_name="$(basename "$r1_chunk")"
    r2_chunk="$sample_chunk_dir/r2/${chunk_name/_R1./_R2.}"
    bam_chunk="$output_dir/bam/${sample_id}.${chunk_name}.bam"
    bwa mem -Y -K 100000000 -t 8 "$reference" "$r1_chunk" "$r2_chunk" \
      | samtools sort -@ 4 -o "$bam_chunk"
    chunk_bams+=("$bam_chunk")
  done

  merged="$output_dir/bam/${sample_id}.merged.bam"
  name_sorted="$output_dir/bam/${sample_id}.name_sorted.bam"
  fixmate="$output_dir/bam/${sample_id}.fixmate.bam"
  coordinate="$output_dir/bam/${sample_id}.coordinate.bam"
  deduplicated="$output_dir/bam/${sample_id}.deduplicated.bam"
  filtered="$output_dir/bam/${sample_id}.chr1.filtered.bam"

  samtools merge -@ 4 -f "$merged" "${chunk_bams[@]}"
  samtools sort -n -@ 4 -o "$name_sorted" "$merged"
  samtools fixmate -m "$name_sorted" "$fixmate"
  samtools sort -@ 4 -o "$coordinate" "$fixmate"
  samtools rmdup -s "$coordinate" "$deduplicated"
  samtools view -@ 4 -b -q 30 -f 3 -F 1796 "$deduplicated" chr1 > "$filtered"
  samtools index "$filtered"
done
