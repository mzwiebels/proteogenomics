#!/usr/bin/Rscript

##########################################################################################
##
## SNV_Signatures.R
##
## Get canonical signatures for samples
##
##########################################################################################

message("\n##Mutation Signature Analysis###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$name)){
  print_help(opt_parser)
  stop("You have to specify a sample name.",call.=FALSE)
} 
if(is.null(opt$species) | !opt$species %in% c("Human","Mouse")){
  print_help(opt_parser)
  stop("You have to specify the species (Human | Mouse).",call.=FALSE)
} 


# main program
if(!require("BSgenome.Hsapiens.UCSC.hg38")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
if(!require("BSgenome.Mmusculus.UCSC.mm10")) BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
if(!require("datasets")) install.packages("datasets")
if(!require("deconstructSigs")) BiocManager::install("deconstructSigs")
if(!require("dplyr")) install.packages("dplyr")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("SomaticCancerAlterations")) BiocManager::install("SomaticCancerAlterations")
if(!require("SomaticSignatures")) BiocManager::install("SomaticSignatures")
if(!require("tidyr")) install.packages("tidyr")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(datasets))
suppressMessages(library(deconstructSigs))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(SomaticCancerAlterations))
suppressMessages(library(SomaticSignatures))
suppressMessages(library(tidyr))

name = opt$name
species = opt$species

file = paste0(name,"/results/Mutect2/",name,".Mutect2.vcf")
sampledf = read.table(file, header=F, sep="\t")
sampledf = sampledf[,c(1,2,4,5)]
colnames(sampledf)=c("chr","pos","ref","alt")
sampledf$Sample = name
sampledf = sampledf[,c("Sample","chr","pos","ref","alt")]

if(species=="Mouse"){
  sampledf = sampledf %>% filter(chr %in% c(1:19,"X","Y"))
  sampledf$chr = paste0("chr",sampledf$chr)
  sigs.input <- mut.to.sigs.input(mut.ref=sampledf,sample.id="Sample",chr="chr",pos="pos",ref="ref",alt="alt",bsg=BSgenome.Mmusculus.UCSC.mm10)
} else if(species=="Human"){
  sampledf = sampledf %>% filter(chr %in% c(1:22,"X","Y"))
  sampledf$chr = paste0("chr",sampledf$chr)
  sigs.input = mut.to.sigs.input(mut.ref=sampledf,sample.id="Sample",chr="chr",pos="pos",ref="ref",alt="alt",bsg=BSgenome.Hsapiens.UCSC.hg38)
}

sample = whichSignatures(tumor.ref=sigs.input,signatures.ref=signatures.nature2013,associated=c("Signature.1A","Signature.2","Signature.3","Signature.4","Signature.5","Signature.6","Signature.7","Signature.8","Signature.9","Signature.10","Signature.11","Signature.12","Signature.13","Signature.14","Signature.15","Signature.16","Signature.17","Signature.18","Signature.19","Signature.20","Signature.21"),sample.id=name,contexts.needed=T,tri.counts.method='default',signature.cutoff=0.2)
pdf(paste0(name,"/results/Mutect2/",name,"_Nature_Pie.pdf",sep=""))
  makePie(sample)
dev.off()
pdf(paste0(name,"/results/Mutect2/",name,"_Nature_Bar.pdf",sep=""))
  plotSignatures(sample)
dev.off()

sample = whichSignatures(tumor.ref=sigs.input,signatures.ref=signatures.cosmic,sample.id=name,contexts.needed=T,tri.counts.method='default',signature.cutoff=0.2)
pdf(paste0(name,"/results/Mutect2/",name,"_Cosmic_Pie.pdf",sep=""))
  makePie(sample)
dev.off()
pdf(paste0(name,"/results/Mutect2/",name,"_Cosmic_Bar.pdf",sep=""))
  plotSignatures(sample)
dev.off()