# Chapter 3.3

required_packages <- c("dplyr", "ggplot2", "pegas", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(ggplot2)
library(pegas)
library(vcfR)

repo_root <- normalizePath(".")
metadata <- read.delim(file.path(repo_root, "data", "metadata", "population_info.tsv"), check.names = FALSE)
locus_manifest <- read.delim(file.path(repo_root, "data", "metadata", "fam72_srgap2_loci.grch38.tsv"), check.names = FALSE)
vcf_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "results", "3.3", "neutrality")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

window_size <- 5000L
step_size <- 500L

window_tajima <- function(vcf, sample_ids, superpopulation, locus) {
  sample_columns <- match(sample_ids, colnames(vcf@gt))
  sample_columns <- sample_columns[!is.na(sample_columns)]
  if (length(sample_columns) == 0) return(data.frame())
  population_vcf <- vcf[, c(1L, sample_columns)]
  positions <- getPOS(population_vcf)
  rows <- list()
  for (start in seq.int(min(positions), max(positions), by = step_size)) {
    end <- start + window_size - 1L
    in_window <- positions >= start & positions <= end
    if (sum(in_window) < 5) next
    dna <- vcfR2DNAbin(population_vcf[in_window, ], extract.indels = FALSE)
    if (length(seg.sites(dna)) < 5) next
    rows[[length(rows) + 1L]] <- data.frame(
      Locus = locus, Superpopulation = superpopulation,
      Window_start = start, Window_end = end, Window_midpoint = (start + end) / 2,
      Segregating_sites = length(seg.sites(dna)), Tajima_D = tajima.test(dna)[[1]]
    )
  }
  bind_rows(rows)
}

all_windows <- bind_rows(lapply(seq_len(nrow(locus_manifest)), function(i) {
  locus <- locus_manifest$locus[i]
  vcf_file <- file.path(vcf_dir, paste0(locus, ".1kg.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing locus VCF: ", vcf_file)
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  bind_rows(lapply(unique(metadata$Superpopulation), function(superpopulation) {
    samples <- metadata$SampleID[metadata$Superpopulation == superpopulation]
    window_tajima(vcf, samples, superpopulation, locus)
  }))
}))

write.table(all_windows, file.path(result_dir, "tajimas_d_fam72_srgap2_5kb_500bp.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot <- ggplot(all_windows, aes(x = Window_midpoint, y = Tajima_D, colour = Superpopulation)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~Locus, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(AFR = "#7F3C8D", AMR = "#11A579", EAS = "#3969AC",
                                 EUR = "#F2B701", SAS = "#E73F74")) +
  labs(x = "GRCh38 chromosome 1 position", y = "Tajima's D", colour = "Superpopulation") +
  theme_bw(base_size = 11) + theme(strip.background = element_blank())
ggsave(file.path(result_dir, "tajimas_d_fam72_srgap2_5kb_500bp.svg"), plot, width = 10, height = 10, units = "in")
