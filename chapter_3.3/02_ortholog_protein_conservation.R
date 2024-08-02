# Chapter 3.3
# The script calculates an entropy-based conservation score and a BLOSUM62-aware 
# score for each column of the vertebrate alignment.
# The former uses BALCONY's position-based sequence weighting, the latter uses
# msa::msaConservationScore with BLOSUM62.

required_packages <- c("BALCONY", "Biostrings", "ggplot2", "msa", "seqinr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(BALCONY)
library(Biostrings)
library(ggplot2)
library(msa)
library(seqinr)
library(tidyr)

repo_root <- normalizePath(".")
alignment_file <- file.path(repo_root, "data", "derived", "3.3", "ortholog_conservation", "fam72a_vertebrates.linsi.fa")
result_dir <- file.path(repo_root, "results", "3.3", "ortholog_conservation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

alignment <- read.alignment(alignment_file, format = "fasta", forceToLower = FALSE)
weights <- get_pos_based_seq_weights(alignment)

entropy_score <- schneider_conservativity(alignment = alignment, weights = weights)

aa_alignment <- readAAMultipleAlignment(alignment_file, format = "fasta")
data(BLOSUM62, package = "msa")
blosum_score <- msaConservationScore(aa_alignment, BLOSUM62, gapVsGap = 0)
blosum_scaled <- blosum_score / max(blosum_score, na.rm = TRUE)

scores <- data.frame(
  Alignment_position = seq_along(entropy_score),
  Shannon_entropy_conservation = entropy_score,
  BLOSUM62_conservation = blosum_score,
  BLOSUM62_conservation_scaled = blosum_scaled
)
write.table(scores, file.path(result_dir, "ortholog_conservation_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- scores[, c("Alignment_position", "Shannon_entropy_conservation", "BLOSUM62_conservation_scaled")]
plot_data <- pivot_longer(plot_data, -Alignment_position,
  names_to = "Metric", values_to = "Conservation_score"
)
plot <- ggplot(plot_data, aes(x = Alignment_position, y = Conservation_score, colour = Metric)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Protein alignment position", y = "Conservation score") +
  theme_bw(base_size = 11)
ggsave(file.path(result_dir, "ortholog_conservation_scores.svg"), plot, width = 10, height = 4, units = "in")
