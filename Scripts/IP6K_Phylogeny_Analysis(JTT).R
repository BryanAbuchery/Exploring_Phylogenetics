#!/usr/bin/env Rscript


############################################################################
# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin17.0 (64-bit)
############################################################################

# This script is based on phylogenetic analysis of IP6K based on amino acid
# sequence fasta file downloaded from NCBI "seqdump"

############################################################################
# Install required packages

# Install the alignment tool (msa). Be sure to check out the documentation on msa;
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("msa")

# Ensure that Tex/LaTex is also installed

#install.packages("pacman")

#install.packages("phangorn")
# read documentation on phangorn use here; https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html

# install.packages("phylotools")
# read documentation on phylotools here; https://cran.r-project.org/web/packages/phylotools/phylotools.pdf

#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("annotate")

#install.packages("rentrez")

#install.packages("seqRFLP")

##########################################################################
# File clean up

# Inspect your sequence file to determine a neat clean up
library("Biostrings")
fastaFile <- readAAStringSet("seqdump.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)


def <- df$seq_name
f1 <- gsub(".*\\[", "", def)
f2 <- gsub("\\].*", "", f1)
df$Organism <- f2
IP6K_for_use_df <- df

#library("seqRFLP") # convert dataframe into fasta file
#df.fasta = dataframe2fas(IP6K_for_use_df, file="IP6K.fasta")

############################################################################
# Load the required packages

library("pacman") # This code chunk shows the loading of required 
# packages. Here, p_load() from pacman, which installs 
# the package if necessary and loads it for use. 
# You can also load installed packages with library() from base R.
pacman::p_load(
  rio,             # import/export
  here,            # relative file paths
  tidyverse,       # general data management and visualization
  ape,             # to import and export phylogenetic files and analyses
  ggtree,          # to visualize phylogenetic files
  treeio,          # to visualize phylogenetic files
  ggnewscale)      # to add additional layers of color schemes
remotes::install_github("YuLab-SMU/tidytree", force = TRUE) #ensures ggtree will work with dplyr
library(msa) # for multiple sequence alignment
library(phangorn)  # phylogeny analyses tool
library(phylotools) # plot and label tree
library(annotate) # Using R environments for annotation
library(rentrez) # Using R to access the NCBI eUtils API

# Construct and write a ∗.txt file with an annotated data frame
# that would include species names and
# sequences with the retrieved sequences, reference for later analyses
write.table(IP6K_for_use_df, file = "IP6K_results_aminoacids.txt", sep = "\t")

# Read this saved file as follows
IP6K_for_use_df <- read.table("IP6K_results_aminoacids.txt", header = TRUE, 
                              sep = "\t", stringsAsFactors = FALSE)

# Construct a file in fasta format that can be used for sequence alignments:
IP6K_definitions <- IP6K_for_use_df$Organism 
extract_names_1_function <- function(x){sub(".∗\\[(.∗)\\].∗", "\\1", x, perl=TRUE)}
taxa_names_seq <- extract_names_1_function(IP6K_definitions)
taxa_names_seq <- gsub(" ", "_", taxa_names_seq)
taxa_fasta_name <- paste0(">", taxa_names_seq)

# Construct fasta entries
taxa_fasta_seq_vector <- paste0(taxa_fasta_name, "\n", IP6K_for_use_df$sequence)

# Remove repeated sequences that have the same accession number:
taxa_fasta_seq_vector_2 <- taxa_fasta_seq_vector[!duplicated(taxa_names_seq)]
taxa_fasta_seq_vector_2 <- taxa_fasta_seq_vector_2[!is.na(taxa_fasta_seq_vector_2)]

# Remove entries with extraneous annotations (e.g., "synthetic_construct", "Chain_A,"):
taxa_fasta_seq_vector_2 <- taxa_fasta_seq_vector_2[!grepl("synthetic_construct|Chain_A", 
                                                          taxa_fasta_seq_vector_2)]

# Write a file in fasta format (see Note 7) with sequences that
# include the annotation for the species of origin and the corresponding accession numbers:
writeLines(taxa_fasta_seq_vector_2, "IP6K_aminoacid.fasta")

###########################################################################
# Alignment

# Uses MSA's Clustal Omega
mySequences <- readAAStringSet("IP6K_aminoacid.fasta") # read fasta file into alignment tool
myFirstAlignment <- msa(mySequences, "ClustalOmega") # align sequences, read documentation for different alignment tools
tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
aldir <- getwd()
workdir <- file.path(aldir,"one_file_per_aln/")
outdir <- file.path(aldir,"alignments")
fasta_file <- paste(outdir,".fasta",sep="")
tex_file <- paste(outdir,".tex",sep="")
msaPrettyPrint(myFirstAlignment, output="tex", showConsensus = "none", askForOverwrite=FALSE, verbose=FALSE,
               file = tex_file, alFile = fasta_file ) #writes the alignment files

###########################################################################
# Maximum Likelihood Tree Inference

# Phylogenetic Inference Using Aligned Sequences in R using distance-based methods
alnfile <- read.aa("alignments.fasta", format="fasta") # read alignment file into phangorn
alnfile_phyDat <- phyDat(alnfile, type = "AA", levels = NULL) # convert to object phylo
aa_dist <- dist.ml(alnfile, model="JTT") # choose suitable model according to your data
alnfile_NJ  <- NJ(aa_dist)

# Save the tree to a file
ape::write.tree(alnfile_NJ, file = "NJ_tree_IP6K_aminoacid.tree")

# Optimize an input tree using a maximum likelihood approach in R.
# Model Testing
mt <- modelTest(alnfile_phyDat, model=c("JTT", "LG", "WAG"), multicore=TRUE)
env <- attr(mt, "env")
fitStart <- eval(get(mt$Model[which.min(mt$BIC)], env), env)
fitStart

# Choose best model according to the "fitStart" (computer chooses)
fitNJ <- pml(alnfile_NJ, alnfile_phyDat, model="JTT", k=4, inv=.2)
fitJTT <- optim.pml(fitNJ, rearrangement = "stochastic", 
                    optInv = TRUE, optGamma = TRUE)
logLik(fitJTT)
bs <- bootstrap.pml(fitJTT, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace = 0))

# Save the tree to a file
ape::write.tree(bs, file = "JTT_tree_IP6K_aminoacid_boot.tree") # writes tree into a newick file

###########################################################################
# Visualizing and Saving Phylogenetic Trees

# The last part of most phylogenetic analyses is the visualization of
# the resulting phylogenies and preparing figures as ∗.jpg or ∗.pdf.

#### Read NJ tree in R and assign them to objects:
IP6K_NJ_tree <- ape::read.tree("NJ_tree_IP6K_aminoacid.tree")

# Open a pdf file and save tree plot:
pdf("IP6K_boostrap_aminoacid_tree.pdf", width=14, height=14)
plot(IP6K_NJ_tree, type = "unrooted", main=" NJ_IP6K aminoacid", cex = 1.2) 
add.scale.bar(cex = 1.0, font = 2, col = "red")
#### Read ML tree in R and assign them to objects:
IP6K_ML_boot_tree <- ape::read.tree("JTT_tree_IP6K_aminoacid_boot.tree")

# best tree with bootstrap values:
phangorn::plotBS(tree = midpoint(IP6K_NJ_tree),
                 BStrees = IP6K_ML_boot_tree,
                 p = 50,
                 type = "fan",
                 bs.col = "blue",
                 frame = "none", cex = 1.2)
title("ML 100 bootstrap replicates based on aligned IP6K amino acid")
add.scale.bar(cex = 0.7, font = 2, col = "red")
dev.off()

###########################################################################
# Save WorkSpace for future use/reference
save.image("Phylogeny_WSP2.RData")