#!/usr/bin/env Rscript


############################################################################
# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin17.0 (64-bit)
############################################################################

# Exploring phylogenetic alignment using DECIPHER

############################################################################
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
#install.packages("seqinr") #to manipulate the sequence file
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
#install.packages("phylotools")
library(phylotools) #Include the tree labels
newtree <- as.phylo(merges)
rename.tips.phylo <- function(tree, names) {
  tree$tip.label <- tree$tip.label <- newtips[match(tree$tip.label,tips)]
  return(tree)
}
tree2plot <- rename.tips.phylo(newtree, newtips)
p <- ggtree(tree2plot, layout = "circular", branch.length = 'none') %<+% sample_data + # %<+% adds dataframe with sample data to tree
  geom_tippoint(
    mapping = aes(color = Organism),          # tip color by Organism. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "Organism",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+
  ggtitle("Phylogenetic tree of Phosphoinositol Pyrophosphates")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 
ggsave("IP6K_tree_circular2.png", plot = p, width = 12, height = 14)
ggsave(plot = p, filename = "Inositol_hexakisphosphate2.pdf",
       width = 12, height = 14, dpi = "retina")
