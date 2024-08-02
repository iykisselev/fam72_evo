#!/usr/bin/env bash
# Chapter 3.3

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_fasta="$repo_root/data/raw/eggnog/fam72a_vertebrates.fa"
output_dir="$repo_root/data/derived/3.3/ortholog_conservation"

mkdir -p "$output_dir"
alignment="$output_dir/fam72a_vertebrates.linsi.fa"

# MAFFT L-INS-i is expressed by --localpair plus --maxiterate 1000.
mafft --localpair --maxiterate 1000 "$input_fasta" > "$alignment"
iqtree2 --sequence "$alignment" --model MFP --bootstrap 1000 --nt AUTO --redo
