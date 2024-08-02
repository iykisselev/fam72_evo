# Chapter 3.3
#
# rehh polarizes alleles with that field, integrates EHH for haplotypes observed at
# least twice, uses an EHH cutoff of 0.01 and 20 kb maximum gap, standardizes iHS
# in 0.01 frequency bins for MAF >= 0.01, and applies Benjamini–Hochberg correction.

required_packages <- c("dplyr", "rehh")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(rehh)

repo_root <- normalizePath(".")
input_dir <- file.path(repo_root, "data", "derived", "3.3", "ihs", "populations")
result_dir <- file.path(repo_root, "results", "3.3", "ihs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

vcf_files <- list.files(input_dir, pattern = "\\.ihs_input\\.vcf\\.gz$", full.names = TRUE)

calculate_population_ihs <- function(vcf_file) {
  population <- sub("\\.ihs_input\\.vcf\\.gz$", "", basename(vcf_file))
  haplotypes <- data2haplohh(
    hap_file = vcf_file,
    polarize_vcf = TRUE,
    vcf_reader = "data.table"
  )
  scan <- scan_hh(
    haplotypes,
    polarized = TRUE,
    limhaplo = 2,
    limehh = 0.01,
    maxgap = 20000,
    interpolate = FALSE,
    discard_integration_at_border = FALSE,
    threads = 1
  )
  ihs <- ihh2ihs(scan, freqbin = 0.01, min_maf = 0.01)
  ihs <- as.data.frame(ihs)
  ihs$Population <- population
  ihs$P_value <- 2 * pnorm(-abs(ihs$IHS))
  ihs$BH_adjusted_P_value <- p.adjust(ihs$P_value, method = "BH")
  ihs$Significant <- ihs$BH_adjusted_P_value < 0.05
  ihs
}

ihs_results <- bind_rows(lapply(vcf_files, calculate_population_ihs))
write.table(ihs_results, file.path(result_dir, "ihs_all_populations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filter(ihs_results, Significant), file.path(result_dir, "ihs_bh_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
