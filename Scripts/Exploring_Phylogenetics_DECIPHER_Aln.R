#!/usr/bin/env Rscript

#ensure that Tex/LaTex is also installed
install.packages("pacman")
library("pacman") # This code chunk shows the loading of required 
# packages. Here, p_load() from pacman, which installs 
# the package if necessary and loads it for use. 
# You can also load installed packages with library() from base R.
pacman::p_load(
  rio,             # import/export
  here,            # relative file paths
  tidyverse,       # general data management and visualization
  ape,             # to import and export phylogenetic files
  ggtree,          # to visualize phylogenetic files
  treeio,          # to visualize phylogenetic files
  ggnewscale)      # to add additional layers of color schemes
remotes::install_github("YuLab-SMU/tidytree", force = TRUE) #ensures ggtree will work with dplyr
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DECIPHER") #for alignment
library(DECIPHER)
AA <- readAAStringSet("concat_IP6K_sequences.fasta", format = "fasta")
uAA <- unique(AA)
index <- match(AA, uAA)
U_AA <- AlignSeqs(uAA)
AA_aligned <- U_AA[index]
writeXStringSet(AA_aligned, file="DECIPHER_Aligned.fasta") #output the alignment file.
install.packages("seqinr") #to manipulate the sequence file
library(seqinr)
Sqr_AA <- read.alignment("DECIPHER_Aligned.fasta", format = "fasta")
distmatrix <- dist.alignment(Sqr_AA, matrix = c("similarity", "identity")) #calculate distance matrix
phylotree <- nj(distmatrix) #create the phylogenetic tree with Neighbour-Joining
ape::write.tree(phylotree, file = "IP6K_tree2.txt")
sample_data <- import("Inositol_hexakisphosphate.csv")
treedf <- as_tibble(phylotree)
tips <- phylotree$tip.label
merges <- merge(treedf, sample_data, by.x = "label", by.y = "Entry", all.x = TRUE, all.y = FALSE)
newtips <- merges$Organism
install.packages("phylotools")
library(phylotools) #Include the tree labels
newtree <- as.phylo(merges)
rename.tips.phylo <- function(tree, names) {
  tree$tip.label <- tree$tip.label <- newtips[match(tree$tip.label,tips)]
  return(tree)
}
tree2plot <- rename.tips.phylo(newtree, newtips)
plots <- ggtree(tree2plot, layout="circular", branch.length = "none", size = .5)
plots2 <- ggtree(tree2plot, branch.length = "none")
ggsave("IP6K_tree_circular2.png", plot = plots)
ggsave("IPK6_tree_linear2.png", plot = plots2)
ggsave(plot = plots, filename = "Inositol_hexakisphosphate(circular)2.pdf",
       width = 10, height = 10, dpi = "retina")
ggsave(plot = plots2, filename = "Inositol_hexakisphosphate(linear)2.pdf",
       width = 10, height = 10, dpi = "retina")



