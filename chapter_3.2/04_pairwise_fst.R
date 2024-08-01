# Chapter 3.2

required_packages <- c("adegenet", "ggplotify", "hierfstat", "pheatmap", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(adegenet)
library(ggplotify)
library(hierfstat)
library(pheatmap)
library(vcfR)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "results", "3.2", "genetic_differentiation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(metadata_file, check.names = FALSE)
if (!all(c("SampleID", "Population", "Superpopulation") %in% names(metadata))) {
  stop("population_info.tsv must contain SampleID, Population, and Superpopulation.")
}

superpopulation_colours <- c(AFR = "#7F3C8D", AMR = "#11A579", EAS = "#3969AC",
                             EUR = "#F2B701", SAS = "#E73F74")

calculate_fst <- function(gene) {
  vcf_file <- file.path(vcf_dir, paste0(gene, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing VCF: ", vcf_file)

  genind <- read.vcfR(vcf_file, verbose = FALSE) %>% vcfR2genind()
  sample_ids <- rownames(genind@tab)
  population <- metadata$Population[match(sample_ids, metadata$SampleID)]
  if (anyNA(population)) stop("A VCF sample is not in population_info.tsv: ", gene)
  genind@pop <- factor(population)

  matrix <- as.matrix(genet.dist(genind, method = "WC84"))
  matrix[matrix < 0] <- 0
  diag(matrix) <- 0
  matrix
}

fst_matrices <- setNames(lapply(paste0("FAM72", LETTERS[1:4]), calculate_fst),
                         paste0("FAM72", LETTERS[1:4]))

saveRDS(fst_matrices, file.path(result_dir, "wc84_pairwise_fst.rds"))

for (gene in names(fst_matrices)) {
  matrix <- fst_matrices[[gene]]
  annotation <- data.frame(
    Superpopulation = metadata$Superpopulation[match(rownames(matrix), metadata$Population)],
    row.names = rownames(matrix)
  )
  long <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
  names(long) <- c("Population_1", "Population_2", "FST")
  long$Gene <- gene
  write.table(long, file.path(result_dir, paste0(tolower(gene), "_wc84_pairwise_fst.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  svg(file.path(result_dir, paste0(tolower(gene), "_wc84_pairwise_fst.svg")), width = 8, height = 7)
  pheatmap(matrix,
    main = paste(gene, "pairwise WC84 FST"),
    annotation_row = annotation,
    annotation_col = annotation,
    annotation_colors = list(Superpopulation = superpopulation_colours),
    border_color = NA
  )
  dev.off()
}
