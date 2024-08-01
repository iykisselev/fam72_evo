# Chapter 3.2
# FAM72 haplotype frequencies by superpopulation

required_packages <- c("dplyr", "ggplot2", "pegas", "stringr", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(ggplot2)
library(pegas)
library(stringr)
library(vcfR)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "results", "3.2", "haplotype_frequencies")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

population_info <- read.delim(metadata_file, check.names = FALSE)
required_columns <- c("SampleID", "Population", "Superpopulation")
if (!all(required_columns %in% names(population_info))) {
  stop("population_info.tsv must contain: ", paste(required_columns, collapse = ", "))
}

gene_names <- paste0("FAM72", LETTERS[1:4])
minimum_haplotype_count <- 50L
palette <- c(AFR = "#7F3C8D", AMR = "#11A579", EAS = "#3969AC",
             EUR = "#F2B701", SAS = "#E73F74")

extract_haplotypes <- function(vcf_file, gene) {
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  dna <- vcfR2DNAbin(vcf, extract.indels = FALSE)
  haps <- haplotype(dna)
  hap_indices <- attr(haps, "index")
  hap_names <- rownames(haps)
  sequence_names <- labels(dna)

  assignments <- bind_rows(lapply(seq_along(hap_indices), function(i) {
    sequence_index <- hap_indices[[i]]
    sequence_id <- sequence_names[sequence_index]
    data.frame(
      Gene = gene,
      Haplotype = hap_names[i],
      Haplotype_sequence = sequence_id,
      SampleID = str_remove(sequence_id, "_[01]$"),
      stringsAsFactors = FALSE
    )
  }))

  assignments <- assignments %>%
    add_count(Haplotype, name = "Haplotype_count") %>%
    filter(Haplotype_count >= minimum_haplotype_count) %>%
    left_join(population_info, by = "SampleID")

  if (anyNA(assignments$Superpopulation)) {
    stop("At least one VCF sample does not occur in population_info.tsv: ", vcf_file)
  }
  assignments
}

assignments <- bind_rows(lapply(gene_names, function(gene) {
  vcf_file <- file.path(vcf_dir, paste0(gene, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing VCF: ", vcf_file)
  extract_haplotypes(vcf_file, gene)
}))

frequency_table <- assignments %>%
  count(Gene, Haplotype, Superpopulation, name = "Haplotype_count") %>%
  group_by(Gene, Haplotype) %>%
  mutate(Total_count = sum(Haplotype_count)) %>%
  ungroup() %>%
  mutate(Haplotype = reorder(Haplotype, Total_count))

write.table(assignments,
  file.path(result_dir, "haplotype_assignments.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(frequency_table,
  file.path(result_dir, "haplotype_frequencies.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot <- ggplot(frequency_table,
  aes(x = Haplotype, y = Haplotype_count, fill = Superpopulation)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Gene, scales = "free_y") +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(x = "Haplotype", y = "Phased-chromosome count", fill = "Superpopulation") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank())

ggsave(file.path(result_dir, "haplotype_frequencies.svg"), plot,
  width = 10, height = 8, units = "in"
)
