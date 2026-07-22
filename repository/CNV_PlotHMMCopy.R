#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotHMMCopy.R
##
## Plot raw data from HMMCopy.
##
##########################################################################################
message("\n##Plot HMMcopy results###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-r", "--repository"),type="character",default="repository",help="path to script repository [default: ./repository/]"),
  make_option(c("-t", "--sequencing_type"),type="character",default=NULL,help="sequencing type (WES | WGS)"),
  make_option(c("-o", "--resolution"),type="integer",default=20000,help="binning resolution [default: 20,000]"),
  make_option(c("-m", "--map_file"),type="character",default=NULL,help="HMMcopy mappability file"),
  make_option(c("-g", "--gc_file"),type="character",default=NULL,help="HMMcopy GC-bias file"),
  make_option(c("-c", "--centromere_file"),type="character",default=NULL,help="centromere region file"),
  make_option(c("-v", "--varregions_file"),type="character",default=NULL,help="variable region file")
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
if(is.null(opt$sequencing_type) | !opt$sequencing_type %in% c("WES","WGS")){
  print_help(opt_parser)
  stop("You have to specify the sequencing type (WES | WGS).",call.=FALSE)
} 

# main program
if(!require("DNAcopy")) install.packages("DNAcopy")
if(!require("GenomeInfoDb")) install.packages("GenomeInfoDb")
if(!require("GenomicRanges")) install.packages("GenomicRanges")
if(!require("HMMcopy")) install.packages("HMMcopy")
if(!require("naturalsort")) install.packages("naturalsort")
suppressMessages(library(DNAcopy))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(GenomicRanges))
suppressMessages(library(HMMcopy))
suppressMessages(library(naturalsort))

name = opt$name
species = opt$species
repository_dir = opt$repository
sequencing_type = opt$sequencing_type
resolution = opt$resolution
map_file = opt$map_file
gc_file = opt$gc_file
centromere_file = opt$centromere_file
varregions_file = opt$varregions_file

# read in wig files and correct for GC and mappability bias
normal = wigsToRangedData(paste0(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig"),gc_file,map_file)
normal$reads = normal$reads+1
normal = as.data.frame(correctReadcount(normal))
normal_copy = GRanges(normal$chr,IRanges(normal$start,normal$end),copy=normal$copy)

tumor = wigsToRangedData(paste0(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig"),gc_file,map_file)
tumor$reads = tumor$reads+1
tumor = as.data.frame(correctReadcount(tumor))
tumor_copy = GRanges(tumor$chr,IRanges(tumor$start,tumor$end),copy=tumor$copy)

# remove regions with increased variability for mice and centromere regions for humams
if (species == "Human"){
	filtering=read.csv(centromere_file,header=F,sep="\t")
	flankLength=5000000
} else if (species == "Mouse"){
	filtering=read.csv(varregions_file,header=F,sep="\t")
	flankLength=0
}

colnames(filtering)[1:3] <- c("space","start","end")
filtering$start = filtering$start-flankLength
filtering$end = filtering$end+flankLength
filtering = GRanges(filtering$space,IRanges(filtering$start,filtering$end))

hits = findOverlaps(query = normal_copy, subject = filtering)
ind = queryHits(hits)
normal_copy = (normal_copy[-ind, ])
message("Removed ",length(ind)," bins near centromeres.")

hits = findOverlaps(query=tumor_copy,subject=filtering)
ind = queryHits(hits)
tumor_copy=(tumor_copy[-ind, ])
message("Removed ", length(ind), " bins near centromeres.")

# computation of the copy number states from the log fold change
somatic_copy = tumor_copy
somatic_copy$copy = tumor_copy$copy-normal_copy$copy
somatic_tab = as.data.frame(somatic_copy)
colnames(somatic_tab) = c("Chrom","Start","End","width","strand","log2Ratio")
somatic_tab = somatic_tab[,c("Chrom","Start","End","log2Ratio")]

write.table(somatic_tab,paste0(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".log2RR.txt"),quote=F,row.names=F,col.names=T,sep='\t')

# segmentation of the CN plot
somatic_CNA = smooth.CNA(CNA(genomdat=somatic_tab$log2Ratio,chrom=somatic_tab$Chrom,maploc=somatic_tab$Start,data.type='logratio'))

if (sequencing_type == "lcWGS"){
	cnv_segments = segment(somatic_CNA,alpha=0.0001,min.width=3,undo.splits='sdundo',undo.SD=1.5,verbose=2)$output
} else if (sequencing_type == "WES" | sequencing_type == "WGS"){
	cnv_segments = segment(somatic_CNA,alpha=0.0001,min.width=5,undo.splits='sdundo',undo.SD=2,verbose=2)$output
}

colnames(cnv_segments) = c("ID","Chrom","Start","End","num.mark","Mean")
cnv_segments = cnv_segments[,c("Chrom","Start","End","Mean")]
cnv_segments = cnv_segments[naturalorder(cnv_segments$Chrom),]

write.table(cnv_segments,paste0(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".segments.txt"),quote=F,row.names=F,col.names=T,sep='\t')

#start plotting
source(paste0(repository_dir,"/all_GeneratePlots.R"))

setwd(paste0(name,"/results/HMMCopy"))
system(paste0("mkdir -p ",name,"_Chromosomes"))
chrom.sizes = DefineChromSizes(species)

if (species=="Human"){
	chromosomes=22
} else if(species=="Mouse"){
	chromosomes=19
}

#define normalization mode, choose from "Mode" or "MAD"
for (y_axis in c("CNV_5","CNV_2")){
	Segments = paste0(name,".HMMCopy.",resolution,".segments.txt")
	Counts = paste0(name,".HMMCopy.",resolution,".log2RR.txt")

	Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="HMMCopy")
	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="HMMCopy")

	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

	for (i in 1:chromosomes){
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
	}
	
	system(paste0("pdfunite ",name,"_Chromosomes/",name,".Chr?.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf ",name,"_Chromosomes/",name,".Chr??.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf ",name,".Chromosomes.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf"))
}