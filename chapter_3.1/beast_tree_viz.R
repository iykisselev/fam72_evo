# Load necessary libraries
library(dplyr)
library(treeio)
library(ggtree)
library(ggplot2)
library(stringr)

# Read BEAST-generated consensus phylogenetic tree
beast <- read.beast("path_to/beast_tree_file")

# Extract sample labels from the tree and restructure the sample data frame
samples <- as_tibble(beast) %>%
  mutate(
    label_extract = str_extract(label, "fam72[a-d]\\.HG[0-9]+"),  # Extract gene and haplotype IDs
    label_parts = str_split_fixed(label_extract, "\\.", 2)  # Split the extracted string at the dot
  ) %>%
  select(-label_extract)  # Remove the temporary column containing haplotype IDs

# Visualize the phylogenetic tree
height <- ggtree(beast) +
  theme_tree2() +
  hexpand(0.2) +
  geom_range("height_0.95_HPD", color='#3D688F', size=3, alpha=0.75) +
  geom_cladelab(
    node = c(31, 38, 44, 46, 29, 4, 1),  # Specify the nodes for labels
    label = c("FAM72A", "FAM72B", "FAM72C", "FAM72D", "P. t. schweinfurthii", 
              "P. t. troglodytes", "Gorilla gorilla"), # Add grouping labels
    align = TRUE,
    fontface = 3,
    offset = -12.5
  )

# Reverse the tree and adjust the x-axis
p <- revts(height) +
  scale_x_continuous(breaks = c(-15:0), labels = abs(-15:0)) +
  theme(plot.margin = unit(c(0, -12, 0, -1), "cm"))

# Save the plot to an SVG file
svg("fam72_neast_tree.svg", width = 8.5, height = 6.5)
print(p) 
dev.off()
