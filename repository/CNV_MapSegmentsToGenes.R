#!/usr/bin/Rscript

##########################################################################################
##
## CNV_MapSegmentsToGenes.R
##
## Takes the segment file from either HMMCopy or Copywriter and maps them to genes.
##
##########################################################################################

message("\n##Filter SNV Output###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-g", "--gencode"),type="character",default=NULL,help="gencode gene annotation file"),
  make_option(c("-m", "--method"),type="character",default=NULL,help="input method (Copywriter | HMMCopy)"),
  make_option(c("-r", "--resolution"),type="integer",default=20000,help="binning resolution [default: 20,000]"),
  make_option(c("-c", "--CGC"),type="character",default=NULL,help="CGC annotation file"),
  make_option(c("-t", "--TruSight"),type="character",default=NULL,help="TruSight annotation file")
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
if(is.null(opt$method) | !opt$method %in% c("Copywriter","HMMCopy")){
  print_help(opt_parser)
  stop("You have to specify the input format (Mutect2 | Strelka).",call.=FALSE)
}


# main program
if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")
if(!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if(!require("tidyr")) install.packages("tidyr")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyr))

name = opt$name
species= opt$species
gencode_file = opt$gencode
method = opt$method
resolution = opt$resolution
CGC = opt$CGC
TruSight = opt$TruSight


genesDT = readRDS(gencode_file)
genesGR = makeGRangesFromDataFrame(genesDT,keep.extra.columns=T)

AnnotateSegment <- function(segDF){
  segDF$chr = paste0("chr",segDF$chr)
  segGR <- makeGRangesFromDataFrame(segDF,keep.extra.columns=T)
  hits <- findOverlaps(segGR,genesGR)
  
  returnDat <- genesDT[subjectHits(hits),c("chr","start","end","geneName","geneID")]
  return(data.frame(returnDat))
}

if (method=="Copywriter"){
	segment = paste0(name,"/results/",method,"/",name,".",method,".segments.Mode.txt")
} else if (method=="HMMCopy") {
	segment = paste0(name,"/results/",method,"/",name,".",method,".",resolution,".segments.txt")
}

cnv = data.frame(Name=NULL,Chrom=NULL,Start=NULL,End=NULL,Mean=NULL,Gene=NULL)
segment = read.delim(segment)

for (i in 1:nrow(segment)){
	temp = NULL
	temp = as.data.frame(segment[i,c("Chrom","Start","End")])
	colnames(temp) = c("chr","start","end")
	results = AnnotateSegment(temp)
	
	if (nrow(results) > 0){
		colnames(results) = c("Chrom","Start","End","Gene","GeneID")
		results$Chrom = gsub("chr","",results$Chrom)
		results$Mean = segment$Mean[i]
		results$Name = name
		results = results[,c("Name","Chrom","Start","End","Mean","Gene","GeneID")]
		cnv = rbind(cnv,results)
	}
}

if (species=="Mouse"){
	segment=makeGRangesFromDataFrame(segment,keep.extra.columns=T)
	ncruc.gr=GRanges(4, IRanges(89311040, 89511040))
	olaps=findOverlaps(segment,ncruc.gr)
	ncruc=data.frame(pintersect(segment[queryHits(olaps)], ncruc.gr[subjectHits(olaps)]))[,c("seqnames","start","end","Mean")]
	ncruc$Gene="Cdkn2_ncruc"
	ncruc$Name=name
	ncruc$GeneID=NA
	colnames(ncruc)=c("Chrom", "Start", "End", "Mean", "Gene", "Name", "GeneID")
	ncruc=ncruc[,c("Name","Chrom", "Start", "End", "Mean", "Gene", "GeneID")]
	cnv=rbind(cnv,ncruc)
}

cnv = tbl_df(cnv) %>% mutate_each(as.character)

cnv = cnv %>% 
  group_by(Gene) %>% 
  arrange(desc(abs(as.numeric(Mean))),.by_group=T) %>% 
  filter(row_number()==1) %>% 
  arrange(as.numeric(Chrom),as.numeric(Start),as.numeric(End))

if (method=="Copywriter"){
	write.table(cnv,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
} else if (method=="HMMCopy"){
	write.table(cnv,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

if(!(is.null(CGC) | CGC=="NULL")){
  CGC=read.delim(CGC,header=T,sep="\t")
  
  cnv_cgc = cnv %>%
    filter(Gene %in% CGC[,1])
  
  if (method=="Copywriter"){
    write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  } else if (method=="HMMCopy") {
    write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  cnv_cgc = cnv_cgc %>%
    filter(abs(as.numeric(Mean)) > 0.75)
  
  if (method=="Copywriter"){
    write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.OnlyImpact.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  } else if (method=="HMMCopy") {
    write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.OnlyImpact.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
}

if (species=="Human" & !(is.null(TruSight) | TruSight=="NULL")){
	TruSight=read.delim(TruSight,header=T,sep="\t")

	cnv_ts = cnv %>%
	filter(Gene %in% TruSight[,1])

	if (method=="Copywriter"){
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	} else if (method=="HMMCopy") {
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	}

	cnv_ts = cnv_ts %>%
	filter(abs(as.numeric(Mean)) > 0.75)

	if (method=="Copywriter"){
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.OnlyImpact.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	} else if (method=="HMMCopy") {
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.OnlyImpact.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	}
}