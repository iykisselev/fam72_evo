# Chapter 3.3 
# Functional classification of gnomAD FAM72 variants

required_packages <- c("BSgenome.Hsapiens.UCSC.hg38", "dplyr", "ggplot2",
                       "GenomicFeatures", "tidyr", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "VariantAnnotation")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R/Bioconductor packages: ", paste(missing_packages, collapse = ", "))

library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)

repo_root <- normalizePath(".")
manifest <- read.delim(file.path(repo_root, "data", "metadata", "fam72.grch38.tsv"), check.names = FALSE)
vcf_dir <- file.path(repo_root, "data", "derived", "3.3", "gnomad")
result_dir <- file.path(repo_root, "results", "3.3", "paralog_conservation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

classify_gene <- function(gene) {
  vcf_file <- file.path(vcf_dir, paste0(gene, ".gnomad_v3.1.2.filtered.vcf.gz"))
  if (!file.exists(vcf_file)) stop("Missing filtered gnomAD VCF: ", vcf_file)

  vcf <- readVcf(vcf_file, genome = "hg38")
  coding <- predictCoding(vcf, TxDb.Hsapiens.UCSC.hg38.knownGene, Hsapiens)
  coding <- as.data.frame(mcols(coding))
  coding$Gene <- gene

  if (nrow(coding) == 0) return(coding)
  coding$Protein_position <- vapply(coding$PROTEINLOC, function(x) {
    if (length(x) == 0 || is.na(x[1])) NA_real_ else as.numeric(x[1])
  }, numeric(1))
  coding
}

annotation <- bind_rows(lapply(manifest$gene, classify_gene))
if (nrow(annotation) == 0) stop("No coding variants were annotated from the filtered gnomAD VCFs.")

write.table(annotation, file.path(result_dir, "gnomad_fam72_coding_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

consequence_counts <- annotation %>% count(Gene, CONSEQUENCE, name = "Variant_count")
write.table(consequence_counts, file.path(result_dir, "gnomad_fam72_consequence_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- annotation %>% filter(!is.na(Protein_position)) %>%
  mutate(Consequence = as.character(CONSEQUENCE))
plot <- ggplot(plot_data, aes(x = Protein_position, y = 0, colour = Consequence)) +
  geom_segment(aes(xend = Protein_position, yend = 1), linewidth = 0.4) +
  geom_point(aes(y = 1), size = 1.6) +
  facet_wrap(~Gene, scales = "free_x") +
  labs(x = "Protein position", y = NULL, colour = "Consequence") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank())
ggsave(file.path(result_dir, "gnomad_fam72_protein_lollipops.svg"), plot, width = 10, height = 5, units = "in")
