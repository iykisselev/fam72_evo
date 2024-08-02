# Chapter 3.3
# BetaScan standardization requires a theta map. This script  calculates Watterson's theta in 
# non-overlapping 0.5 Mb chr1 windows for each 1000 Genomes population and writes the three-column format required by BetaScan.py:
# zero-based window start, exclusive window end, and theta.

required_packages <- c("dplyr", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(vcfR)

repo_root <- normalizePath(".")
vcf_file <- file.path(repo_root, "data", "derived", "3.3", "ihs", "chr1.1kg.ancestral_allele_snps.vcf.gz")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
output_dir <- file.path(repo_root, "data", "derived", "3.3", "betascan", "theta_maps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(metadata_file, check.names = FALSE)
if (!file.exists(vcf_file)) stop("Missing ancestral-allele VCF: ", vcf_file)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
positions <- getPOS(vcf)
window_width <- 500000L

allele_count <- function(genotype_vector) {
  genotype_vector <- genotype_vector[!is.na(genotype_vector) & genotype_vector != "./."]
  alleles <- unlist(strsplit(gsub("[|/]", "", genotype_vector), ""))
  length(unique(alleles))
}

theta_for_population <- function(population) {
  sample_ids <- metadata$SampleID[metadata$Population == population]
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)[, sample_ids, drop = FALSE]
  windows <- data.frame(
    start = seq.int(1L, max(positions), by = window_width),
    end = pmin(seq.int(1L, max(positions), by = window_width) + window_width, max(positions) + 1L)
  )
  values <- lapply(seq_len(nrow(windows)), function(i) {
    in_window <- positions >= windows$start[i] & positions < windows$end[i]
    gt_window <- gt[in_window, , drop = FALSE]
    segregating_sites <- sum(apply(gt_window, 1, function(x) allele_count(x) > 1L))
    non_missing_haplotypes <- median(rowSums(!is.na(gt_window) & gt_window != "./.") * 2, na.rm = TRUE)
    if (!is.finite(non_missing_haplotypes) || non_missing_haplotypes < 2) {
      return(data.frame(start = windows$start[i] - 1L, end = windows$end[i] - 1L, theta = NA_real_))
    }
    a1 <- sum(1 / seq_len(non_missing_haplotypes - 1))
    # BetaScan expects theta per base pair. Divide Watterson's window estimate by
    # the actual window width before writing the theta map.
    theta <- if (segregating_sites > 0 && a1 > 0) {
      (segregating_sites / a1) / (windows$end[i] - windows$start[i])
    } else {
      NA_real_
    }
    data.frame(start = windows$start[i] - 1L, end = windows$end[i] - 1L, theta = theta)
  })
  bind_rows(values)
}

for (population in unique(metadata$Population)) {
  theta_map <- theta_for_population(population) %>% filter(!is.na(theta))
  write.table(theta_map, file.path(output_dir, paste0(population, ".theta_map.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
