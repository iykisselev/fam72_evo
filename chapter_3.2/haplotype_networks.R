# Execute the following code within an RStudio session.
# Load necessary libraries
library(vcfR)
library(dplyr)
library(pegas)
library(stringr)

# Function to generate a network plot from VCF data
network_plt <- function(vcf_filename, popdata){
  # Read VCF data and convert to DNAbin format
  dnabin <- read.vcfR(vcf_filename) %>% vcfR2DNAbin() 
  
  # Extract haplotypes with a minimum frequency threshold
  h <- haplotype(dnabin) %>% subset(minfreq = 15)
  
  # Create a haplotype network
  nt <- haploNet(h)
  
  # Retrieve labels and sizes for plotting
  nt_labs <- attr(nt, "labels")
  sz <- summary(h)[nt_labs]
  
  # Match sample IDs after stripping indices and retrieve corresponding superpopulation
  pop <- popdata$Superpopulation[match(str_remove(labels(dnabin), "_[0-1]"), 
                                       popdata$SampleID)]
  
  # Calculate haplotype frequencies based on superpopulation
  P <- haploFreq(dnabin, fac = pop, haplo = h)[nt_labs, ]
  
  # Set colors for the superpopulations
  setHaploNetOptions(pie.colors.function = c("#7F3C8D", "#11A579", "#3969AC", 
                                             "#F2B701", "#E73F74"))
  
  # Plot the network
  plt <- plot(nt, size = sz, pie = P, labels = FALSE, legend = TRUE,
              show.mutation = 0, scale.ratio = 200)
}

# Read population data
popdata <- read.delim("20130606_g1k_3202_samples_ped_population.txt", sep = " ")
fam72_grch38 <- read.table("fam72.grch38.tsv", header = FALSE)

# Generate network plots for each FAM72 gene
lapply(fam72_grch38$V4[1:4], function(gene) {
network_plt(paste0(gene, ".1KG.vcf.gz"), popdata)
})
