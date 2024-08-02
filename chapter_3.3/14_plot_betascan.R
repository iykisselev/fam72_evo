# Chapter 3.3 
# Identify and plot extreme standardized beta(1) values

required_packages <- c("dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(ggplot2)

repo_root <- normalizePath(".")
score_dir <- file.path(repo_root, "results", "3.3", "betascan")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
manifest_file <- file.path(repo_root, "data", "metadata", "fam72.grch38.tsv")

score_files <- list.files(score_dir, pattern = "\\.betascores\\.txt$", full.names = TRUE)
metadata <- read.delim(metadata_file, check.names = FALSE)
manifest <- read.delim(manifest_file, check.names = FALSE)

read_scores <- function(score_file) {
  population <- sub("\\.betascores\\.txt$", "", basename(score_file))
  scores <- read.delim(score_file, check.names = FALSE)
  scores$Population <- population
  scores
}

scores <- bind_rows(lapply(score_files, read_scores)) %>%
  group_by(Population) %>%
  mutate(Extreme = Beta1_std >= quantile(Beta1_std, 0.98, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(distinct(metadata, Population, Superpopulation), by = "Population")

locus_scores <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  filter(scores, Position >= manifest$start[i] & Position <= manifest$end[i]) %>%
    mutate(Gene = manifest$gene[i])
}))

write.table(filter(locus_scores, Extreme), file.path(score_dir, "betascan_beta1_top_2_percent_fam72.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot <- ggplot(locus_scores, aes(x = Position, y = Beta1_std)) +
  geom_point(colour = "grey75", size = 0.5) +
  geom_point(data = filter(locus_scores, Extreme), aes(colour = Superpopulation), size = 1.1) +
  facet_wrap(~Gene, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(AFR = "#7F3C8D", AMR = "#11A579", EAS = "#3969AC",
                                 EUR = "#F2B701", SAS = "#E73F74")) +
  labs(x = "GRCh38 chromosome 1 position", y = "Standardized beta(1)", colour = "Superpopulation") +
  theme_bw(base_size = 11) + theme(strip.background = element_blank())
ggsave(file.path(score_dir, "betascan_beta1_fam72.svg"), plot, width = 9, height = 9, units = "in")
