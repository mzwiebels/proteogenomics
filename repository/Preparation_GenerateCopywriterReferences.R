#!/usr/bin/Rscript

##########################################################################################
##
## Preparation_GenerateCopywriterReferences.R
##
## Creates Copywriter reference files for 20kb windows.
##
##########################################################################################
library("CopywriteR")
args <- commandArgs(TRUE)

species = args[1]
version = args[2]

ref.genome = ifelse(species=="Mouse","mm10","hg38")

for (resolution in c(10000,20000,50000,100000)) {
	preCopywriteR(output.folder = paste0("ref/",version),bin.size = resolution,ref.genome = ref.genome)
}