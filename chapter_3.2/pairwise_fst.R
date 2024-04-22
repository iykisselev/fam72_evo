# This code reads genomic data from FAM72A-D VCF files, calculates pairwise Fst, and visualizes the results as heatmaps. It processes genomic regions specified in the "fam72.grch38.tsv" file, stores the results in a list, saves the Fst matrices, and exports the combined heatmap plots to an SVG file in pop_diff filder.
library(vcfR)
library(adegenet)
library(pheatmap)
library(ggpubr)
library(ggplotify)
library(hierfstat)
library(RColorBrewer)

# Read gene coordinate data
fam72_grch38 <- read.table("fam72.grch38.tsv", header = FALSE)

# Function to calculate pairwise FST for specified genomic regions in a VCF file
wc_pairwise_fst <- function(vcf_file, start, end){
  fam72_vcfr <- read.vcfR(vcf_file)
  popdata <- read.delim("20130606_g1k_3202_samples_ped_population.txt", sep=" ")
  fam72_vcfr_clean <- fam72_vcfr[getPOS(fam72_vcfr) >= start & getPOS(fam72_vcfr) <= end,]
  fam72_genind <- vcfR2genind(fam72_vcfr_clean)
  fam72_genind@pop <- as.factor(popdata$Population[match(rownames(fam72_genind@tab), 
                                                         popdata$SampleID)])
  fam72_pairwise_fst <- genet.dist(fam72_genind, method = "WC84")
  fam72_pairwise_fst[fam72_pairwise_fst < 0] <- 0  # Ensure no negative FST values
  return(fam72_pairwise_fst)
}

# Function to create a heatmap from FST matrix
fst_pheatmap <- function(pairwise_fst, title){
  data <- as.matrix(pairwise_fst)
  popdata <- read.delim("20130606_g1k_3202_samples_ped_population.txt", sep=" ") # Load population data
  popanot <- data.frame(Superpopulation = popdata$Superpopulation[match(rownames(data), popdata$Population)], 
                        row.names=rownames(data))
  mycolors <- c("AFR" = "#7F3C8D", "AMR" = "#11A579", "EAS" = "#3969AC", "EUR" = "#F2B701", "SAS" = "#E73F74")
  p <- pheatmap(data, main = title, annotation_col = popanot, annotation_row = popanot, 
                annotation_colors = list(Superpopulation = mycolors),
                breaks = seq(0, 0.25, by = 0.0025))
  return(as.ggplot(p))
}

# Initialize list to store FST matrices
fst_matrices <- vector(mode = "list", length = nrow(fam72_grch38))

# Calculate pairwise FST for each region using data from fam72_grch38 and generate heatmaps
results <- lapply(1:nrow(fam72_grch38), function(i) {
  vcf_filename <- paste0(fam72_grch38$V4[i], ".1KG.vcf.gz")
  fst <- wc_pairwise_fst(vcf_filename, fam72_grch38$V2[i], fam72_grch38$V3[i])
  heatmap_plot <- fst_pheatmap(fst, paste("FAM72", LETTERS[i], sep = ""))
  list(fst_matrix = fst, heatmap = heatmap_plot)
})

# Populate fst_matrices with results from lapply and gather all plots
plots <- vector("list", length = length(results))
names(fst_matrices) <- paste("FAM72", LETTERS[1:length(results)], sep = "")

for (i in seq_along(results)) {
  fst_matrices[[i]] <- results[[i]]$fst_matrix
  plots[[i]] <- results[[i]]$heatmap
}

# Arrange plots in a grid layout
combined_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

# Create a directory for the results of population differentiation analysis
dir.create("pop_diff")

# Save Fst matrices for the nMDS 
saveRDS(fst_matrices, "./pop_diff/wc_pairwise_fst.Rdata")

# Save the plots
svg("./pop_diff/fam72_beast_tree.svg", width = 20, height = 13)
print(combined_plot)
dev.off()
