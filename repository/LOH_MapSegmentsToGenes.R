#!/usr/bin/Rscript

##########################################################################################
##
## LOH_MapSegmentsToGenes.R
##
## Takes the segment file from Titan and maps it to genes.
##
##########################################################################################

message("\n##Map Titan Segments to Genes###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-g", "--gencode_file"),type="character",default=NULL,help="Gencode annotation file"),
  make_option(c("-c", "--CGC"),type="character",default=NULL,help="CGC annotation file")
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
if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")
if(!require("GenomicRanges")) install.packages("GenomicRanges")
if(!require("tidyr")) install.packages("tidyr")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyr))

name = opt$name
species= opt$species
gencode_file = opt$gencode_file
CGC = opt$CGC


genesDT = readRDS(gencode_file)
genesGR <- makeGRangesFromDataFrame(genesDT, keep.extra.columns = T)

AnnotateSegment <- function(segDF){
  segDF[,"chr"]=paste0("chr",segDF[,"chr"])
  segGR <- makeGRangesFromDataFrame(segDF, keep.extra.columns = T)
  hits <- findOverlaps(segGR, genesGR)
  
  #returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneID)] # only geneID
  returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneName, geneID)] # more stuff
  
  return(data.frame(returnDat))
}

segments = paste(name,"/results/Titan/run_ploidy2/",name,"_cluster01.segs.txt",sep="")
loh = data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL,TITAN=NULL,Gene=NULL)
segments = tbl_df(read.delim(segments))

segments = segments %>%
  filter(TITAN_call %in% c("ALOH","NLOH","DLOH","HOMD")) %>%
  filter(Length.snp. >= 10) %>%
  filter((End_Position.bp.-Start_Position.bp.)/Length.snp. < 1000000)

if (nrow(segments) > 0){
	for (i in 1:nrow(segments)){
		temp = as.data.frame(segments[i,c("Chromosome","Start_Position.bp.","End_Position.bp.")])
		colnames(temp) = c("chr","start","end")
		
		results = AnnotateSegment(temp)
		if (nrow(results) > 0){
			results = results[,c("chr", "start", "end", "geneName","geneID")]
			colnames(results) = c("Chrom", "Start", "End", "Gene","GeneID")
			results[,"Chrom"] = gsub("chr","",results[,"Chrom"])
			results$TITAN = as.data.frame(segments)[i,"TITAN_call"]
			results$Name = name
			results = results[,c("Name", "Chrom", "Start", "End", "TITAN","Gene","GeneID")]
			loh = rbind(loh,results)
		}
	}

	loh = tbl_df(loh) %>% mutate_each(as.character)
	loh = loh %>% 
	  #group_by(Gene) %>% 
	  #arrange(TITAN,.by_group=T) %>% 
	  #filter(row_number()==1) %>% 
	  arrange(as.numeric(Chrom),as.numeric(Start),as.numeric(End))
} else{
	loh[1,"Name"] = name
	loh[1,"Chrom"] = 1
	loh[1,"Start"] = 1
	loh[1,"End"] = 1
	loh[1,"TITAN"] = "NLOH"
	loh[1,"Gene"] = "EMPTY"
	loh[1,"GeneID"] = "EMPTY"
}

write.table(loh,paste(name,"/results/LOH/",name,".LOH.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

if(!is.null(CGC)){
  CGC = read.delim(CGC,header=T,sep="\t")
  loh_cgc = loh %>% filter(Gene %in% as.character(CGC[,1]))
  
  write.table(loh_cgc,paste(name,"/results/LOH/",name,".LOH.genes.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}