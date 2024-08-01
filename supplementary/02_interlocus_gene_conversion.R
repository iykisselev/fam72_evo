# This script filters a published interlocus-gene-conversion table to FAM72 events,
# maps donor and acceptor intervals to the FAM72 paralog coordinates and plots them

required_packages <- c("circlize", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(circlize)
library(dplyr)
library(readr)

repo_root <- normalizePath(".")
igc_file <- file.path(repo_root, "data", "raw", "supplementary", "interlocus_gene_conversion.csv")
manifest_file <- file.path(repo_root, "data", "metadata", "fam72.grch38.tsv")
result_dir <- file.path(repo_root, "results", "supplementary", "interlocus_gene_conversion")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(igc_file)) stop("Missing interlocus-gene-conversion table: ", igc_file)

manifest <- read.delim(manifest_file, check.names = FALSE)
events <- read_csv(igc_file, show_col_types = FALSE) %>%
  filter(grepl("FAM72", IGC.genes)) %>%
  select(chrm.acceptor, start.acceptor, end.acceptor,
         chrm.donor, start.donor, end.donor, Seen.in.n.samples)

assign_paralog <- function(start, end) {
  hit <- which(start <= manifest$end & end >= manifest$start)
  if (length(hit) == 0) NA_character_ else manifest$gene[hit[1]]
}

events$Acceptor_gene <- mapply(assign_paralog, events$start.acceptor, events$end.acceptor)
events$Donor_gene <- mapply(assign_paralog, events$start.donor, events$end.donor)
events <- filter(events, !is.na(Acceptor_gene), !is.na(Donor_gene))
write.table(events, file.path(result_dir, "fam72_interlocus_gene_conversion_events.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sectors <- manifest %>% transmute(chr = gene, start = 0, end = end - start + 1)
donor <- events %>% transmute(chr = Donor_gene, start = 1, end = pmax(2, end.donor - start.donor + 1))
acceptor <- events %>% transmute(chr = Acceptor_gene, start = 1, end = pmax(2, end.acceptor - start.acceptor + 1))

svg(file.path(result_dir, "fam72_interlocus_gene_conversion.svg"), width = 8, height = 8)
circos.clear()
circos.par(gap.degree = 8)
circos.genomicInitialize(sectors, sector.width = rep(1, nrow(sectors)))
circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 0.5, CELL_META$sector.index, facing = "bending.inside", cex = 0.8)
})
circos.genomicLink(donor, acceptor, col = add_transparency("#3D688F", 0.65), border = NA)
dev.off()
