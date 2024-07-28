# Chapter 3.1
# SHAPEIT v2 cannot phase multiple records with the same position. Records created
# by splitting multiallelic sites are assigned consecutive temporary positions and
# the original coordinate is retained in INFO/ORIG_POS. 03_restore_...R reverses
# this change after phasing.

required_packages <- c("VariantAnnotation", "Rsamtools")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(VariantAnnotation)

repo_root <- normalizePath(".")
vcf_dir <- file.path(repo_root, "data", "derived", "3.1", "chimpanzee", "vcf")
vcf_files <- list.files(vcf_dir, pattern = "\\.filtered\\.vcf\\.gz$", full.names = TRUE)
if (length(vcf_files) == 0) stop("No filtered chimpanzee VCFs found in ", vcf_dir)

for (vcf_file in vcf_files) {
  vcf <- readVcf(vcf_file)
  original_positions <- start(rowRanges(vcf))
  occurrence <- ave(original_positions, original_positions, FUN = seq_along)
  info(vcf)$ORIG_POS <- as.integer(original_positions)
  start(rowRanges(vcf)) <- original_positions + occurrence - 1L
  end(rowRanges(vcf)) <- start(rowRanges(vcf)) + nchar(as.character(ref(vcf))) - 1L

  output_file <- sub("\\.filtered\\.vcf\\.gz$", ".shifted.vcf.gz", vcf_file)
  writeVcf(vcf, output_file, index = TRUE)
}
