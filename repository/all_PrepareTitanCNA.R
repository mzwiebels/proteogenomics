#!/usr/bin/Rscript

##########################################################################################
##
## all_RunTitanCNA.R
##
## Prepares all files for running Titan.
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
  make_option(c("-e", "--exons_file"),type="character",default=NULL,help="exon region file")
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
if(!require("HMMcopy")) install.packages("HMMcopy")
if(!require("GenomeInfoDb")) install.packages("GenomeInfoDb")
if(!require("GenomicRanges")) install.packages("GenomicRanges")
if(!require("TitanCNA")) install.packages("TitanCNA")
suppressMessages(library(DNAcopy))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(GenomicRanges))
suppressMessages(library(HMMcopy))

name = opt$name
species = opt$species
repository_dir = opt$repository
sequencing_type = opt$sequencing_type
resolution = opt$resolution
map_file = opt$map_file
gc_file = opt$gc_file
exons_file = opt$exons_file


correctReadDepth = function(tumWig,normWig,gcWig,mapWig,genomeStyle="NCBI",targetedSequence=NULL){
  message("Reading GC and mappability files")
  gc = wigToGRanges(gcWig)
  map = wigToGRanges(mapWig)
  
  message("Loading tumor file: ", tumWig)
  tumor_reads = wigToGRanges(tumWig)
  message("Loading normal file: ", normWig)
  normal_reads = wigToGRanges(normWig)
  
  seqlevelsStyle(gc) = genomeStyle
  seqlevelsStyle(map) = genomeStyle
  seqlevelsStyle(tumor_reads) = genomeStyle
  seqlevelsStyle(normal_reads) = genomeStyle
  
  gc = gc[seqnames(gc) %in% seqnames(tumor_reads)]
  map = map[seqnames(map) %in% seqnames(tumor_reads)]
  samplesize = 50000
  
  if(!is.null(targetedSequence)){
    message("Analyzing targeted regions...")
    targetIR = GRanges(ranges=IRanges(start=targetedSequence[,2],end=targetedSequence[,3]),seqnames=targetedSequence[,1])
    names(targetIR) = setGenomeStyle(seqlevels(targetIR),genomeStyle)
    
    hits = findOverlaps(query=tumor_reads,subject=targetIR)
    keepInd = unique(queryHits(hits))
    
    tumor_reads = tumor_reads[keepInd,]
    normal_reads = normal_reads[keepInd,]
    gc = gc[keepInd,]
    map = map[keepInd,]
    samplesize = min(ceiling(nrow(tumor_reads)*0.1),samplesize)
  }
  tumor_reads$gc = gc$value
  tumor_reads$map = map$value
  colnames(values(tumor_reads)) = c("reads","gc","map")
  normal_reads$gc = gc$value
  normal_reads$map = map$value
  colnames(values(normal_reads)) = c("reads","gc","map")
  
  message("Correcting Tumor")
  tumor_copy = correctReadcount(tumor_reads,samplesize=samplesize)
  
  message("Correcting Normal")
  normal_copy = correctReadcount(normal_reads,samplesize=samplesize)
  
  message("Normalizing Tumor by Normal")
  tumor_copy$copy = tumor_copy$copy-normal_copy$copy
  rm(normal_copy)
  
  temp = cbind(chr=as.character(seqnames(tumor_copy)),start=start(tumor_copy),end=end(tumor_copy),logR=tumor_copy$copy)
  temp = as.data.frame(temp,stringsAsFactors=F)
  
  mode(temp$start) = "numeric"
  mode(temp$end) = "numeric"
  mode(temp$logR) = "numeric"
  
  return(temp)
}

# read in wig files and correct for GC and mappability bias
tumWig = paste0(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig")
normWig = paste0(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig")
variants = paste0(name,"/results/LOH/",name,".VariantsForLOH.txt")

if(sequencing_type == "WGS"){
	cnData = correctReadDepth(tumWig,normWig,gc_file,map_file,genomeStyle="NCBI") 
} else if(sequencing_type == "WES"){
	exons_file = read.delim(exons_file)
	exons_file = exons_file[,1:3]
	colnames(exons_file) = c("chr","start","end")
	cnData = correctReadDepth(tumWig,normWig,gc_file,map_file,genomeStyle="NCBI",exons_file) 
}

cnData[is.na(cnData$logR),"logR"] = 0
write.table(cnData,paste0(name,"/results/Titan/",name,".cnFile.txt"),sep="\t",quote=F,row.names=F,col.names=T)

variants = read.delim(variants)
variants = variants[variants$Normal_Freq <= 0.7 & variants$Normal_Freq >= 0.3,]
data = variants[,c("Chrom","Pos","Ref","Tumor_RefCount","Alt","Tumor_AltCount")]

write.table(data,paste(name,"/results/Titan/",name,".hetFile.txt",sep=""),sep="\t", quote=F, row.names=F,col.names=T)