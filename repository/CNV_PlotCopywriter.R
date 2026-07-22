#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotCopywriter.R
##
## Plot raw data from Copywriter.
##
##########################################################################################
message("\n##Plot CopywriteR results###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-r", "--repository"),type="character",default="repository",help="path to script repository [default: ./repository/]"),
  make_option(c("-m", "--normalization"),type="character",default="Mode",help="normalization strategy (Mode | MAD) [default: Mode]")
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
if(is.null(opt$normalization) | !opt$normalization %in% c("Mode","MAD")){
  print_help(opt_parser)
  stop("You have to specify the species (Human | Mouse).",call.=FALSE)
} 

# main program
name = opt$name
species = opt$species
repository_dir = opt$repository
normalization = opt$normalization

source(paste0(repository_dir,"/all_GeneratePlots.R"))

setwd(paste0(name,"/results/Copywriter"))
system(paste0("mkdir -p ",name,"_Chromosomes"))

chrom.sizes = DefineChromSizes(species)

if (species=="Human"){
	chromosomes=22
} else if (species=="Mouse"){
	chromosomes=19
}

#define normalization mode, choose from "Mode" or "MAD"
for (y_axis in c("CNV_5","CNV_2")){
	Segments = paste(name,".Copywriter.segments.",normalization,".txt",sep="")
	Counts = paste(name,".Copywriter.log2RR.",normalization,".txt",sep="")

	Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="Copywriter")
	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="Copywriter")

	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="Copywriter",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

	for (i in 1:chromosomes){
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="Copywriter",normalization=normalization,y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7,outformat="pdf")
	}

	system(paste0("pdfunite ",name,"_Chromosomes/",name,".Chr?.CNV.Copywriter.",normalization,".",gsub("CNV_","",y_axis),".pdf ",name,"_Chromosomes/",name,".Chr??.CNV.Copywriter.",normalization,".",gsub("CNV_","",y_axis),".pdf ",name,".Chromosomes.CNV.Copywriter.",normalization,".",gsub("CNV_","",y_axis),".pdf"))
}