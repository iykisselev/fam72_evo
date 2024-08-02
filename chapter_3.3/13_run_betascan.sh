#!/usr/bin/env bash
# Calculate standardized beta(1) in 2 kb windows

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
betascan_program="$repo_root/tools/BetaScan/BetaScan.py"
input_dir="$repo_root/data/derived/3.3/betascan/betascan_input"
theta_dir="$repo_root/data/derived/3.3/betascan/theta_maps"
result_dir="$repo_root/results/3.3/betascan"

mkdir -p "$result_dir"

for input_file in "$input_dir"/*.betascan.txt.gz; do
  [[ -s "$input_file" ]] || continue
  population="$(basename "$input_file" .betascan.txt.gz)"
  theta_map="$theta_dir/${population}.theta_map.tsv"
  [[ -s "$theta_map" ]] || { echo "Missing theta map: $theta_map" >&2; exit 1; }
  python3 "$betascan_program" \
    -i "$input_file" \
    -o "$result_dir/${population}.betascores.txt" \
    -w 2000 \
    -std \
    -theta_map "$theta_map"
done
