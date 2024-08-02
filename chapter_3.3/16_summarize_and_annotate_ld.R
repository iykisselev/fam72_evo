# Chapter 3.3 
# The script identifies top-two-percent LD pairs within each population and locus, extracts the first position of each pair, counts the
# populations supporting each position, retains the gene-wise top ten percent of
# positions, bins them into 100 bp intervals, and intersects them with the Ensembl
# Regulatory Build release 110 GFF3

required_packages <- c("dplyr", "GenomicRanges", "rtracklayer")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R/Bioconductor packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(GenomicRanges)
library(rtracklayer)

repo_root <- normalizePath(".")
pair_dir <- file.path(repo_root, "data", "derived", "3.3", "ld_pairs")
regulatory_gff <- file.path(repo_root, "data", "raw", "ensembl", "regulatory_features_grch38.gff3.gz")
result_dir <- file.path(repo_root, "results", "3.3", "ld")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

pair_files <- list.files(pair_dir, pattern = "_ld_pairs\\.tsv$", full.names = TRUE)
if (length(pair_files) == 0) stop("No LD pair tables found. Run 15_compute_ld.R first.")
pairs <- bind_rows(lapply(pair_files, read.delim, check.names = FALSE))

top_pairs <- pairs %>%
  group_by(Locus, Population) %>%
  mutate(Top_2_percent = LD_r_squared >= quantile(LD_r_squared, 0.98, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Top_2_percent)

position_support <- top_pairs %>%
  distinct(Locus, Population, SNP_1) %>%
  count(Locus, SNP_1, name = "Supporting_populations") %>%
  group_by(Locus) %>%
  mutate(Top_10_percent = Supporting_populations >= quantile(Supporting_populations, 0.90, na.rm = TRUE)) %>%
  ungroup()

candidate_intervals <- position_support %>%
  filter(Top_10_percent) %>%
  mutate(Bin_start = ((SNP_1 - 1L) %/% 100L) * 100L + 1L,
         Bin_end = Bin_start + 99L) %>%
  distinct(Locus, Bin_start, Bin_end, .keep_all = TRUE)

write.table(top_pairs, file.path(result_dir, "ld_top_2_percent_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(position_support, file.path(result_dir, "ld_position_population_support.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidate_intervals, file.path(result_dir, "ld_top_10_percent_100bp_intervals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

regulatory_features <- import(regulatory_gff)
candidate_ranges <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = candidate_intervals$Bin_start, end = candidate_intervals$Bin_end),
  Locus = candidate_intervals$Locus,
  Supporting_populations = candidate_intervals$Supporting_populations
)
overlaps <- findOverlaps(candidate_ranges, regulatory_features, ignore.strand = TRUE)

annotation <- cbind(
  as.data.frame(candidate_ranges[queryHits(overlaps)]),
  as.data.frame(regulatory_features[subjectHits(overlaps)])
)
write.table(annotation, file.path(result_dir, "ld_candidate_regulatory_overlaps.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
