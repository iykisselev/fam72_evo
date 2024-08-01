required_packages <- c("dplyr", "ggplot2", "ggtranscript", "rtracklayer", "stringr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R/Bioconductor packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(ggplot2)
library(ggtranscript)
library(rtracklayer)
library(stringr)
library(tidyr)

repo_root <- normalizePath(".")
input_dir <- file.path(repo_root, "data", "raw", "supplementary", "gtex_long_read")
result_dir <- file.path(repo_root, "results", "supplementary", "long_read_transcripts")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

tissue_file <- file.path(input_dir, "tissue_annot.csv")
tpm_file <- file.path(input_dir, "quantification_flair_filter.tpm.txt")
gtf_file <- file.path(input_dir, "flair_filter_transcripts.gtf")

target_genes <- c("ENSG00000196550", "ENSG00000188610", "ENSG00000263513", "ENSG00000215784")
tissue <- read.csv(tissue_file, stringsAsFactors = FALSE) %>%
  select(sample_id, tissue) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "."),
         Organ = str_split_fixed(tissue, " - ", 2)[, 1]) %>%
  filter(!Organ %in% c("", "K562"))
tpm <- read.delim(tpm_file, check.names = FALSE)
gtf <- import(gtf_file)

target_transcripts <- unique(gtf$transcript_id[gtf$gene_id %in% target_genes])
expression <- tpm %>%
  filter(transcript %in% target_transcripts) %>%
  pivot_longer(-transcript, names_to = "sample_id", values_to = "TPM") %>%
  left_join(tissue, by = "sample_id") %>%
  filter(!is.na(Organ))

write.table(expression, file.path(result_dir, "long_read_fam72_expression.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

expression_plot <- ggplot(expression, aes(x = Organ, y = TPM)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.45, size = 0.7) +
  facet_wrap(~transcript, scales = "free_y", ncol = 1) +
  labs(x = "Organ", y = "TPM") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank())
ggsave(file.path(result_dir, "long_read_fam72_expression.svg"), expression_plot,
       width = 10, height = 10, units = "in")

exons <- as.data.frame(gtf) %>%
  filter(gene_id %in% target_genes, type == "exon") %>%
  select(seqnames, start, end, strand, gene_id, transcript_id)
structure_plot <- ggplot(exons, aes(xstart = start, xend = end, y = transcript_id)) +
  geom_range() +
  geom_intron(data = to_intron(exons, "transcript_id"), aes(strand = strand)) +
  facet_wrap(~gene_id, scales = "free_x", ncol = 1) +
  labs(x = "Genomic position", y = "Transcript") +
  theme_bw(base_size = 11) + theme(strip.background = element_blank())
ggsave(file.path(result_dir, "long_read_fam72_transcript_structures.svg"), structure_plot,
       width = 10, height = 9, units = "in")
