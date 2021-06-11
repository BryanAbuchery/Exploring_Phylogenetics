#!/usr/bin/env Rscript

#ensure that Tex/LaTex is also installed
#install.packages("pacman")
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
library(DECIPHER)
AA <- readAAStringSet("concat_IP6K_sequences.fasta", format = "fasta")
uAA <- unique(AA)
index <- match(AA, uAA)
U_AA <- AlignSeqs(uAA)
writeXStringSet(AA_aligned, file="DECIPHER_Aligned.fasta")
library(seqinr)
Sqr_AA <- read.alignment("DECIPHER_Aligned.fasta", format = "fasta")
distmatrix <- dist.alignment(Sqr_AA, matrix = c("similarity", "identity"))
phylotree <- nj(distmatrix)
ape::write.tree(phylotree, file = "IP6K_tree2.txt")
sample_data <- import("Inositol_hexakisphosphate.csv")
treedf <- as_tibble(phylotree)
tips <- phylotree$tip.label
merges <- merge(treedf, sample_data, by.x = "label", by.y = "Entry", all.x = TRUE, all.y = FALSE)
newtips <- merges$Organism
library(phylotools)
newtree <- as.phylo(merges)
newtree$tip.label <- newtips
plots <- ggtree(newtree, branch.length = "none", mrsd = "2013-01-01", size = .5)
ggsave("IP6K_tree_linear1.png", plot = plots)
ggsave(plot = plots, filename = "Inositol_hexakisphosphate2.pdf",
       width = 10, height = 10, dpi = "retina")




