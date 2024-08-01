# Chapter 3.2 — PERMANOVA of FST distances
# Two models are fit for each FST matrix: five 1000 Genomes superpopulations and 
# an Out-of-Africa grouping (AFR versus non-AFR). Each model uses 1e6 permutations

required_packages <- "vegan"
if (!requireNamespace(required_packages, quietly = TRUE)) {
  stop("Install required R package: vegan")
}
library(vegan)

repo_root <- normalizePath(".")
metadata_file <- file.path(repo_root, "data", "metadata", "population_info.tsv")
fst_file <- file.path(repo_root, "results", "3.2", "genetic_differentiation", "wc84_pairwise_fst.rds")
result_dir <- file.path(repo_root, "results", "3.2", "genetic_differentiation")

metadata <- read.delim(metadata_file, check.names = FALSE)
fst_matrices <- readRDS(fst_file)

extract_model_row <- function(fit, gene, model) {
  table <- as.data.frame(fit)
  data.frame(
    Gene = gene,
    Model = model,
    R2 = table$R2[1],
    F = table$F[1],
    P_value = table$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )
}

model_results <- do.call(rbind, lapply(names(fst_matrices), function(gene) {
  fst <- fst_matrices[[gene]]
  populations <- rownames(fst)
  design <- metadata[match(populations, metadata$Population), ]
  if (anyNA(design$Superpopulation)) stop("FST populations are not all represented in population_info.tsv.")
  design$OoA <- ifelse(design$Superpopulation == "AFR", "AFR", "OoA")

  superpopulation_fit <- adonis2(as.dist(fst) ~ Superpopulation,
                                  data = design, permutations = 1000000)
  out_of_africa_fit <- adonis2(as.dist(fst) ~ OoA,
                               data = design, permutations = 1000000)
  rbind(
    extract_model_row(superpopulation_fit, gene, "Superpopulation"),
    extract_model_row(out_of_africa_fit, gene, "Out_of_Africa")
  )
}))

best_models <- do.call(rbind, lapply(split(model_results, model_results$Gene), function(gene_results) {
  both_significant <- all(gene_results$P_value < 0.05)
  selected <- if (both_significant) gene_results$Model[which.max(gene_results$R2)] else NA_character_
  data.frame(Gene = gene_results$Gene[1], Both_models_significant = both_significant,
             Selected_model = selected, stringsAsFactors = FALSE)
}))

write.table(model_results, file.path(result_dir, "permanova_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(best_models, file.path(result_dir, "permanova_model_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
