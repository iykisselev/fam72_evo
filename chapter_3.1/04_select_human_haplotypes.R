# Chapter 3.1

required_packages <- c("dplyr", "pegas", "stringr", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(pegas)
library(stringr)
library(vcfR)

repo_root <- normalizePath(".")
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "data", "derived", "3.1", "beast")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

select_haplotypes <- function(gene) {
  vcf_file <- file.path(vcf_dir, paste0(gene, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing VCF: ", vcf_file)
  dna <- read.vcfR(vcf_file, verbose = FALSE) %>% vcfR2DNAbin(extract.indels = FALSE)
  haps <- haplotype(dna)
  indices <- attr(haps, "index")
  labels <- labels(dna)
  haplotype_names <- rownames(haps)

  bind_rows(lapply(seq_along(indices), function(i) {
    members <- labels[indices[[i]]]
    representative <- members[1]
    data.frame(
      Gene = gene,
      Haplotype = haplotype_names[i],
      Haplotype_count = length(members),
      SampleID = str_remove(representative, "_[01]$"),
      Haplotype_copy = as.integer(str_extract(representative, "[01]$")) + 1L,
      stringsAsFactors = FALSE
    )
  })) %>% filter(Haplotype_count > 100)
}

selected <- bind_rows(lapply(paste0("FAM72", LETTERS[1:4]), select_haplotypes))
write.table(selected, file.path(result_dir, "human_haplotypes_above_100.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
