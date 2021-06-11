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
#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("msa")
library(msa)
mySequences <- readAAStringSet("concat_IP6K_sequences.fasta")
myFirstAlignment <- msa(mySequences)
tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
aldir <- getwd()
workdir <- file.path(aldir,"one_file_per_aln/")
outdir <- file.path(aldir,"alignments")
fasta_file <- paste(outdir,".fasta",sep="")
tex_file <- paste(outdir,".tex",sep="")
msaPrettyPrint(myFirstAlignment, output="tex", showConsensus = "none", askForOverwrite=FALSE, verbose=FALSE,
               file = tex_file, alFile = fasta_file )
treeAln <- msaConvert(myFirstAlignment, type="seqinr::alignment")
#install.packages("seqinr")
library(seqinr)
distmatrix <- dist.alignment(treeAln, "identity")
phylotree <- nj(distmatrix)
ape::write.tree(phylotree, file='IP6K_tree.txt')
tips <- phylotree$tip.label
sample_data <- import("Inositol_hexakisphosphate.csv")
treedtf <- as_tibble(phylotree)
mergeddf <- merge(treedtf, sample_data, by.x = "label", by.y = "Entry", all.x =TRUE, all.y = FALSE)
newtips <- mergeddf$Organism
library(phylotools)
new_tree <- as.phylo(mergeddf)
new_tree$tip.label <- newtips
plots <- ggtree(new_tree, branch.length = "none", mrsd = "2013-01-01", size = .5) %<+% mergeddf +
  geom_tiplab(aes(label = Organism), linesize = .5)
ggsave("IP6K_tree_linear.png", plot = plots)
ggsave(plot = plots, filename = "Inositol_hexakisphosphate.pdf",
       width = 10, height = 10, dpi = "retina")
