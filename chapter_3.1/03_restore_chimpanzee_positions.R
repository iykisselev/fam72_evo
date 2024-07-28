# Chapter 3.1
# Restore INFO/ORIG_POS values written by 02_shift_...R.

required_packages <- "VariantAnnotation"
if (!requireNamespace(required_packages, quietly = TRUE)) stop("Install required R package: VariantAnnotation")
library(VariantAnnotation)

repo_root <- normalizePath(".")
phased_dir <- file.path(repo_root, "data", "derived", "3.1", "chimpanzee", "phased")
vcf_files <- list.files(phased_dir, pattern = "\\.shifted\\.phased\\.vcf\\.gz$", full.names = TRUE)
if (length(vcf_files) == 0) stop("No shifted phased chimpanzee VCFs found in ", phased_dir)

for (vcf_file in vcf_files) {
  vcf <- readVcf(vcf_file)
  original_positions <- info(vcf)$ORIG_POS
  if (is.null(original_positions)) stop("INFO/ORIG_POS is absent from ", vcf_file)
  start(rowRanges(vcf)) <- as.integer(original_positions)
  end(rowRanges(vcf)) <- start(rowRanges(vcf)) + nchar(as.character(ref(vcf))) - 1L
  output_file <- sub("\\.shifted\\.phased\\.vcf\\.gz$", ".restored.vcf.gz", vcf_file)
  writeVcf(vcf, output_file, index = TRUE)
}
