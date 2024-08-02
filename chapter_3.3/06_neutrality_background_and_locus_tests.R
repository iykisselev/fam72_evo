# Chapter 3.3 
# Background and locus neutrality statistics

required_packages <- c("dplyr", "pegas", "PopGenome", "vcfR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(pegas)
library(PopGenome)
library(vcfR)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
segment_dir <- file.path(repo_root, "data", "derived", "3.3", "neutrality", "background_vcfs")
locus_dir <- file.path(repo_root, "data", "derived", "3.2", "loci")
result_dir <- file.path(repo_root, "results", "3.3", "neutrality")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

window_width <- 17387L
metadata <- read.delim(metadata_file, check.names = FALSE)

superpopulation_samples <- split(metadata$SampleID, metadata$Superpopulation)

make_windows <- function(start, end, width = window_width) {
  starts <- seq.int(start, end, by = width)
  data.frame(Window_start = starts, Window_end = pmin(starts + width - 1L, end))
}

is_polymorphic <- function(genotypes) {
  alleles <- unlist(strsplit(gsub("[|/]", "", genotypes[!is.na(genotypes) & genotypes != "./."], perl = TRUE), ""))
  length(unique(alleles)) > 1L
}

diversity_for_window <- function(vcf, samples) {
  sample_columns <- match(samples, colnames(vcf@gt))
  sample_columns <- sample_columns[!is.na(sample_columns)]
  if (length(sample_columns) == 0) return(c(Nucleotide_diversity = NA_real_, Haplotype_diversity = NA_real_))
  subset_vcf <- vcf[, c(1L, sample_columns)]
  dna <- vcfR2DNAbin(subset_vcf, extract.indels = FALSE)
  c(Nucleotide_diversity = nuc.div(dna), Haplotype_diversity = hap.div(dna, variance = FALSE, method = "Nei"))
}

popgenome_statistics <- function(vcf_file, population_samples, windows, start, end) {
  # PopGenome receives named superpopulation sample sets, applies the required non-overlapping window transform, and returns the standard neutrality matrix.
  genome <- readVCF(vcf_file, tid = "chr1", frompos = start, topos = end, numcols = 10000)
  genome <- set.populations(genome, population_samples, diploid = TRUE)
  genome <- sliding.window.transform(genome, width = window_width, jump = window_width, type = 2)
  genome <- neutrality.stats(genome, detail = TRUE)
  statistics <- get.neutrality(genome)
  bind_rows(lapply(seq_along(statistics), function(i) {
    population_statistics <- as.data.frame(statistics[[i]])
    if (nrow(population_statistics) != nrow(windows)) {
      stop("PopGenome returned a different number of windows than expected for ", vcf_file)
    }
    population_statistics$Superpopulation <- names(population_samples)[i]
    population_statistics$Window_start <- windows$Window_start
    population_statistics$Window_end <- windows$Window_end
    population_statistics
  }))
}

background_rows <- list()
segment_files <- list.files(segment_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
if (length(segment_files) == 0) stop("No background VCF segments found. Run 05_prepare_neutrality_input.sh first.")

for (segment_file in segment_files) {
  vcf <- read.vcfR(segment_file, verbose = FALSE)
  positions <- getPOS(vcf)
  if (length(positions) == 0) next
  windows <- make_windows(min(positions), max(positions))

  # PopGenome supplies Tajima.D, Fay.Wu.H, and Zeng.E for each superpopulation.
  pg_stats <- popgenome_statistics(segment_file, superpopulation_samples, windows, min(positions), max(positions))
  write.table(pg_stats, file.path(result_dir, paste0(basename(segment_file), ".popgenome_neutrality.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  for (superpopulation in names(superpopulation_samples)) {
    sample_ids <- superpopulation_samples[[superpopulation]]
    for (window_index in seq_len(nrow(windows))) {
      window <- windows[window_index, ]
      row_index <- positions >= window$Window_start & positions <= window$Window_end
      if (sum(row_index) == 0) next
      gt <- extract.gt(vcf[row_index, ], element = "GT", as.numeric = FALSE)
      sample_gt <- gt[, intersect(sample_ids, colnames(gt)), drop = FALSE]
      segregating_sites <- sum(apply(sample_gt, 1, is_polymorphic))
      if (segregating_sites < 5) next
      diversity <- diversity_for_window(vcf[row_index, ], sample_ids)
      pg_row <- pg_stats %>% filter(
        Superpopulation == superpopulation,
        Window_start == window$Window_start,
        Window_end == window$Window_end
      )
      background_rows[[length(background_rows) + 1L]] <- cbind(
        data.frame(Segment = basename(segment_file), Segregating_sites_checked = segregating_sites),
        pg_row,
        data.frame(
          Nucleotide_diversity = diversity[["Nucleotide_diversity"]],
          Haplotype_diversity = diversity[["Haplotype_diversity"]]
        )
      )
    }
  }
}

background <- bind_rows(background_rows)
write.table(background, file.path(result_dir, "neutrality_background_diversity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Empirical five-percent cut-offs are calculated separately for every statistic and  superpopulation. 
metrics <- c("Tajima.D", "Fay.Wu.H", "Zeng.E", "Nucleotide_diversity", "Haplotype_diversity")
thresholds <- bind_rows(lapply(metrics, function(metric) {
  background %>% group_by(Superpopulation) %>% summarise(
    Metric = metric,
    Lower_5_percent = quantile(.data[[metric]], 0.05, na.rm = TRUE),
    Upper_95_percent = quantile(.data[[metric]], 0.95, na.rm = TRUE),
    .groups = "drop"
  )
}))
write.table(thresholds, file.path(result_dir, "neutrality_empirical_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

locus_files <- list.files(locus_dir, pattern = "^FAM72[A-D]\\.1kg\\.vcf\\.gz$", full.names = TRUE)
write.table(data.frame(Locus_VCF = locus_files), file.path(result_dir, "neutrality_locus_input_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
