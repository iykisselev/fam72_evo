# Chapter 3.3
# For every locus and 1000 Genomes population, the script converts the biallelic SNP VCF 
# to DNAbin and calculates all pairwise LD correlations with pegas::LDscan.

required_packages <- c("pegas", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(pegas)
library(vcfR)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
locus_manifest_file <- file.path(repo_root, "data", "metadata", "fam72_srgap2_loci.grch38.tsv")
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
output_dir <- file.path(repo_root, "data", "derived", "3.3", "ld_pairs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(metadata_file, check.names = FALSE)
locus_manifest <- read.delim(locus_manifest_file, check.names = FALSE)

calculate_population_ld <- function(vcf, population, locus) {
  samples <- metadata$SampleID[metadata$Population == population]
  sample_columns <- match(samples, colnames(vcf@gt))
  sample_columns <- sample_columns[!is.na(sample_columns)]
  if (length(sample_columns) < 2) return(data.frame())

  population_vcf <- vcf[, c(1L, sample_columns)]
  dna <- vcfR2DNAbin(population_vcf, extract.indels = FALSE)
  segregating <- seg.sites(dna)
  if (length(segregating) < 2) return(data.frame())

  ld_matrix <- as.matrix(LDscan(dna))
  positions <- getPOS(population_vcf)[segregating]
  if (nrow(ld_matrix) != length(positions)) stop("LD matrix and segregating-site positions do not agree.")
  pairs <- which(lower.tri(ld_matrix), arr.ind = TRUE)
  r <- ld_matrix[pairs]
  data.frame(
    Locus = locus,
    Population = population,
    SNP_1 = positions[pairs[, 1]],
    SNP_2 = positions[pairs[, 2]],
    LD_r = r,
    LD_r_squared = r^2
  )
}

for (i in seq_len(nrow(locus_manifest))) {
  locus <- locus_manifest$locus[i]
  vcf_file <- file.path(vcf_dir, paste0(locus, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing locus VCF: ", vcf_file)
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  for (population in unique(metadata$Population)) {
    pairs <- calculate_population_ld(vcf, population, locus)
    if (nrow(pairs) == 0) next
    write.table(pairs, file.path(output_dir, paste0(locus, "_", population, "_ld_pairs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
