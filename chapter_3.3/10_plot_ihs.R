# Chapter 3.3

required_packages <- c("dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(ggplot2)

repo_root <- normalizePath(".")
ihs_file <- file.path(repo_root, "results", "3.3", "ihs", "ihs_all_populations.tsv")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
manifest_file <- file.path(repo_root, "data", "metadata", "fam72.grch38.tsv")
result_dir <- file.path(repo_root, "results", "3.3", "ihs")

ihs <- read.delim(ihs_file, check.names = FALSE)
metadata <- read.delim(metadata_file, check.names = FALSE)
manifest <- read.delim(manifest_file, check.names = FALSE)

ihs <- ihs %>% left_join(distinct(metadata, Population, Superpopulation), by = "Population")

locus_scores <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  filter(ihs, POSITION >= manifest$start[i] & POSITION <= manifest$end[i]) %>%
    mutate(Gene = manifest$gene[i])
}))

write.table(filter(locus_scores, Significant), file.path(result_dir, "ihs_bh_significant_fam72.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot <- ggplot(locus_scores, aes(x = POSITION, y = IHS)) +
  geom_point(colour = "grey75", size = 0.6) +
  geom_point(data = filter(locus_scores, Significant), aes(colour = Superpopulation), size = 1.2) +
  facet_wrap(~Gene, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(AFR = "#7F3C8D", AMR = "#11A579", EAS = "#3969AC",
                                 EUR = "#F2B701", SAS = "#E73F74")) +
  labs(x = "GRCh38 chromosome 1 position", y = "Standardized iHS", colour = "Superpopulation") +
  theme_bw(base_size = 11) + theme(strip.background = element_blank())
ggsave(file.path(result_dir, "ihs_fam72_paralogs.svg"), plot, width = 9, height = 9, units = "in")
