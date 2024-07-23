# The script calculates and plots frequencies of FAM72A-D haplotypes 
# Load necessary libraries
library(vcfR)
library(dplyr)
library(ggplot2)
library(rcartocolor)
library(ggpubr)
library(pegas)
library(stringr)

# Read population data and genetic coordinates for FAM72 genes
popdata <- read.delim("20130606_g1k_3202_samples_ped_population.txt", sep = " ")
fam72_grch38 <- read.table("fam72.grch38.tsv", header = FALSE)

# Function to extract haplotypes with minimal frequency of 'minfreq' and sort
extract_haplotypes <- function(vcf_file, minfreq) {
  v <- read.vcfR(vcf_file)
  v_dnabin <- vcfR2DNAbin(v, extract.indels = FALSE)
  haps_50 <- haplotype(v_dnabin) %>% subset(minfreq = minfreq) %>% sort()
  
  df_list <- lapply(seq_along(attr(haps_50, "index")), function(i) {
    data.frame(
      Haplotype = rep(attr(haps_50, "dimnames")[[1]][i], 
                      length(attr(haps_50, "index")[[i]])),
      Sample = attr(v_dnabin, "dimnames")[[1]][attr(haps_50, "index")[[i]]]
    )
  })
  
  df <- do.call(rbind, df_list)
  return(df)
}

# Function to annotate haplotypes with populations and superpopulations
process_haplotypes <- function(df, popdata) {
  df <- df %>% mutate(Sample = str_remove_all(Sample, "_.+")) %>%
    left_join(popdata[,c(2,6,7)], by = c("Sample" = "SampleID")) %>%
    group_by(Haplotype, Superpopulation) %>%
    summarise(n = n(), .groups = 'drop')
  
  df$Haplotype <- factor(df$Haplotype, levels = 
                           (df %>% group_by(Haplotype) %>% summarise(f = sum(n)) %>% 
                              arrange(desc(f)))$Haplotype)
  
  return(df)
}

#Function to plot frequencies of haplotypes by superpopulations
plot_haplotypes <- function(df, gene) {
  plt <- ggplot(data = df, aes(y = n, x = Haplotype, fill = Superpopulation)) + 
    geom_bar(position = "stack", stat = "identity") + theme_bw() + 
    scale_fill_manual(values = carto_pal(6, "Bold")[-6]) + ggtitle(gene) + theme(
      legend.text = element_text(size=14), legend.title = element_text(size=14),
      axis.title = element_text(size=12)) + ylab("Haplotype count")
  
  return(plt)
}

# Wrapper function
haplotype_bar <- function(vcf_file, minfreq, gene, popdata) {
  df <- extract_haplotypes(vcf_file, minfreq)
  df_hist <- process_haplotypes(df, popdata)
  plt <- plot_haplotypes(df_hist, gene)
  return(plt)
}

# Generate bar plots for each FAM72 gene based on its genomic region
plots <- lapply(1:nrow(fam72_grch38), function(i) {
  gene_info <- fam72_grch38[i, ]
  vcf_filename <- paste0(gene_info$V4, ".1KG.vcf.gz") 
  haplotype_bar(vcf_filename, minfreq = 50, toupper(gene_info$V4), popdata)
})

# Arrange the haplotype plots in a grid
haplotypes_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

# Save the plot
svg("./pop_diff/hap_freq.svg", width = 11.3, height = 10)
print(haplotypes_plot)
dev.off()
