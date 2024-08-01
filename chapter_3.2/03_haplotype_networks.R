# Chapter 3.2
# FAM72 haplotype networks are inferred with pegas from indel-free phased 1000 Genomes locus VCFs

required_packages <- c("pegas", "stringr", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(pegas)
library(stringr)
library(vcfR)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "results", "3.2", "haplotype_networks")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

population_info <- read.delim(metadata_file, check.names = FALSE)
if (!all(c("SampleID", "Superpopulation") %in% names(population_info))) {
  stop("population_info.tsv must contain SampleID and Superpopulation.")
}

network_palette <- c("#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74")
setHaploNetOptions(pie.colors.function = network_palette)

plot_network <- function(gene) {
  vcf_file <- file.path(vcf_dir, paste0(gene, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing VCF: ", vcf_file)

  dna <- read.vcfR(vcf_file, verbose = FALSE) %>% vcfR2DNAbin(extract.indels = FALSE)
  haplotypes <- haplotype(dna) %>% subset(minfreq = 15)
  network <- haploNet(haplotypes)
  labels_in_network <- attr(network, "labels")
  sizes <- summary(haplotypes)[labels_in_network]

  sample_ids <- str_remove(labels(dna), "_[01]$")
  superpopulations <- population_info$Superpopulation[
    match(sample_ids, population_info$SampleID)
  ]
  if (anyNA(superpopulations)) stop("A VCF sample is not in population_info.tsv: ", gene)

  pie_table <- haploFreq(dna, fac = superpopulations, haplo = haplotypes)[labels_in_network, , drop = FALSE]

  svg(file.path(result_dir, paste0(tolower(gene), "_network.svg")), width = 8, height = 7)
  plot(network, size = sizes, pie = pie_table, labels = FALSE, legend = TRUE,
       show.mutation = 0, scale.ratio = 200, main = gene)
  dev.off()
}

invisible(lapply(paste0("FAM72", LETTERS[1:4]), plot_network))
