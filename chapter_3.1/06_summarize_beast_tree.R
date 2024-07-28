# Chapter 3.1

required_packages <- c("dplyr", "ggplot2", "ggtree", "treeio")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)

repo_root <- normalizePath(".")
tree_file <- file.path(repo_root, "results", "3.1", "beast", "selected_clock.mcc.tree")
result_dir <- file.path(repo_root, "results", "3.1", "beast")
if (!file.exists(tree_file)) stop("Missing TreeAnnotator MCC tree: ", tree_file)

tree <- read.beast(tree_file)
node_table <- as_tibble(tree)
write.table(node_table, file.path(result_dir, "selected_clock_mcc_node_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot <- ggtree(tree) +
  theme_tree2() +
  hexpand(0.15) +
  geom_tiplab(size = 3, align = TRUE) +
  labs(x = "Time before present (million years)")

# height_0.95_HPD is written by BEAST/TreeAnnotator when present in the MCC tree.
if ("height_0.95_HPD" %in% names(node_table)) {
  plot <- plot + geom_range("height_0.95_HPD", colour = "#3D688F", linewidth = 1.1, alpha = 0.7)
}

ggsave(file.path(result_dir, "selected_clock_mcc_tree.svg"), plot,
       width = 10, height = 7, units = "in")
