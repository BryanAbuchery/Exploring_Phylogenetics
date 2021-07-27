#!/usr/bin/env Rscript

#ensure that Tex/LaTex is also installed
#install.packages("pacman", repos = "http://cran.us.r-project.org")
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
#BiocManager::install("msa") #for alignment
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
               file = tex_file, alFile = fasta_file ) #writes the alignment files
treeAln <- msaConvert(myFirstAlignment, type="seqinr::alignment") #convert alignment into a file recognized by seqinr
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
#install.packages("phylotools") #plot and label tree
library(phylotools)
new_tree <- as.phylo(mergeddf)
rename.tips.phylo <- function(tree, names) {
  tree$tip.label <- tree$tip.label <- newtips[match(tree$tip.label,tips)]
  return(tree)
}
tree2plot <- rename.tips.phylo(new_tree, newtips)
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
ggsave("IP6K_tree_circular.png", plot = p, width = 12, height = 14)
ggsave(plot = p, filename = "Inositol_hexakisphosphate.pdf",
       width = 12, height = 14, dpi = "retina")
