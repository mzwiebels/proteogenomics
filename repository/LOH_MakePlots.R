#!/usr/bin/Rscript

##########################################################################################
##
## LOH_MakePlots.R
##
## Plot raw data for LOH vizualisation.
##
##########################################################################################

message("\n##Plot LOH###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-r", "--repository"),type="character",default="repository",help="path to scipt repository  [default: ./repository/]")
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
name = opt$name
species = opt$species
repository_dir = opt$repository

source(paste0(repository_dir,"/all_GeneratePlots.R"))


setwd(paste0(name,"/results/LOH"))
system(paste0("mkdir -p ",name,"_Chromosomes"))

chrom.sizes = DefineChromSizes(species)

if (species=="Human"){
	chromosomes = 22
} else if (species=="Mouse"){
	chromosomes = 19
}

data = paste0(name,".VariantsForLOH.txt")
LOHDat = ProcessCountData(data,chrom.sizes,"LOH")
plotGlobalRatioProfile(cn=LOHDat[[1]],ChromBorders=LOHDat[[2]],cnSeg="",samplename=name,method="LOH",toolname="LOH",normalization="",y_axis="LOH",Transparency=70,Cex=0.3,outformat="pdf")

for(i in 1:chromosomes){
  plotChromosomalRatioProfile(cn=LOHDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH",toolname="LOH",SliceStart="",SliceStop="",Transparency=70,Cex=0.7,outformat="pdf")
}
system(paste0("pdfunite ",name,"_Chromosomes/",name,".Chr?.LOH.LOH.pdf ",name,"_Chromosomes/",name,".Chr??.LOH.LOH.pdf ",name,".Chromosomes.LOH.LOH.pdf"))

LOHDat = ProcessCountData(data,chrom.sizes,"LOH_raw")
plotGlobalRatioProfile(cn=LOHDat[[1]],ChromBorders=LOHDat[[2]],cnSeg="",samplename=name,method="LOH_raw",toolname="LOH_raw",normalization="",y_axis="LOH_raw",Transparency=70,Cex=0.3,outformat="pdf")

for(i in 1:chromosomes){
  plotChromosomalRatioProfile(cn=LOHDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH_raw",toolname="LOH_raw",SliceStart="",SliceStop="",Transparency=70,Cex=0.7, outformat="pdf")
}
system(paste0("pdfunite ",name,"_Chromosomes/",name,".Chr?.LOH_raw.LOH_raw.pdf ",name,"_Chromosomes/",name,".Chr??.LOH_raw.LOH_raw.pdf ",name,".Chromosomes.LOH_raw.LOH_raw.pdf"))


data = paste0(name,".VariantsForLOHGermline.txt")
LOH_GermlineDat = ProcessCountData(data,chrom.sizes,"LOH_Germline")
plotGlobalRatioProfile(cn=LOH_GermlineDat[[1]],ChromBorders=LOH_GermlineDat[[2]],cnSeg="",samplename=name,method="LOH_Germline",toolname="LOH_Germline",normalization="",y_axis="LOH_Germline",Transparency=70, Cex=0.3,outformat="pdf")

for(i in 1:chromosomes){
  plotChromosomalRatioProfile(cn=LOH_GermlineDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH_Germline",toolname="LOH_Germline",SliceStart="",SliceStop="",Transparency=70,Cex=0.7,outformat="pdf")
}
system(paste0("pdfunite ",name,"_Chromosomes/",name,".Chr?.LOH_Germline.LOH_Germline.pdf ",name,"_Chromosomes/",name,".Chr??.LOH_Germline.LOH_Germline.pdf ",name,".Chromosomes.LOH_Germline.LOH_Germline.pdf"))