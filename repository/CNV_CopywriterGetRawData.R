#!/usr/bin/Rscript

##########################################################################################
##
## CNV_CopywriterGetRawData.R
##
## Extract datapoints and segments from the Rdata object provided by CopywriteR.
##
##########################################################################################
message("\n##Extract segments from CopywriteR output###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-r", "--runmode"),type="character",default=NULL,help="runmode (MS | SS)"),
  make_option(c("-t", "--type"),type="character",default=NULL,help="sample type (Tumor | Normal)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$name)){
  print_help(opt_parser)
  stop("You have to specify a sample name.",call.=FALSE)
} 
if(is.null(opt$runmode) | !opt$runmode %in% c("MS","SS")){
  print_help(opt_parser)
  stop("You have to specify the species (Human | Mouse).",call.=FALSE)
} 


# main program
name = opt$name
runmode = opt$runmode
type = opt$type

# Descend into Copywriter results folder
setwd(paste0(name,"/results/Copywriter"))

#load raw segment data from file, change around colunms and export them
load("CNAprofiles/segment.Rdata")
segmentData = segment.CNA.object$output

if (runmode == "MS") {
	Selection = unique(grep("Normal",grep("Tumor",segmentData$ID,value=T),value=T))
} else if (runmode == "SS") {
	Selection = paste0("log2.",gsub("-",".",name),".",type,".bam.vs.none")
} 

segmentData = segmentData[segmentData$ID==Selection,c("chrom","loc.start","loc.end","seg.mean")]
colnames(segmentData)=c("Chrom","Start","End","Mean")
segmentData$Start = floor(segmentData$Start)

write.table(segmentData,file=paste0(name,".Copywriter.segments.MAD.txt"),quote=F,sep="\t",row.names=F,col.names=T)

#calculate the mode of segment means and subtract it from the MAD
d = density(segmentData$Mean[segmentData$Chrom<=22])
Shift = d$x[which.max(d$y)]
segmentData$Mean = segmentData$Mean - Shift

write.table(segmentData,file=paste0(name,".Copywriter.segments.Mode.txt"),quote=F,sep="\t",row.names=F,col.names=T)

#extract raw count read counts
logReadCounts = read.table("CNAprofiles/log2_read_counts.igv",header=T,sep="\t")

if (runmode == "MS") {
  TumorLogReadCounts = grep("Tumor",colnames(logReadCounts),value=T)
  NormalLogReadCounts = grep("Normal",colnames(logReadCounts),value=T)
  logReadCounts$Copy = logReadCounts[,TumorLogReadCounts]-logReadCounts[,NormalLogReadCounts]
  logReadCountsMode = logReadCounts
  
  logReadCountsMode[,TumorLogReadCounts] = logReadCountsMode[,TumorLogReadCounts]-Shift
  logReadCountsMode$Copy = logReadCountsMode[,TumorLogReadCounts]-logReadCountsMode[,NormalLogReadCounts]
  logReadCountsMode = logReadCountsMode[,c("Chromosome","Start","End","Copy")]
  colnames(logReadCountsMode) = c("Chrom","Start","End","log2Ratio")
  write.table(logReadCountsMode,file=paste0(name,".Copywriter.log2RR.Mode.txt"),quote=F,sep="\t",row.names=F,col.names=T)
  
  logReadCounts = logReadCounts[,c("Chromosome","Start","End","Copy")]
  colnames(logReadCounts) = c("Chrom","Start","End","log2Ratio")
  write.table(logReadCounts,file=paste0(name,".Copywriter.log2RR.MAD.txt"),quote=F,sep="\t",row.names=F,col.names=T)
} else if (runmode == "SS") {
  Selection = paste0("log2.",gsub("-",".",name),".",type,".bam")
  
  logReadCounts$Copy = logReadCounts[,Selection]
  logReadCountsMode = logReadCounts
  logReadCountsMode[,Selection] = logReadCountsMode[,Selection]-Shift
  logReadCountsMode = logReadCountsMode[,c("Chromosome","Start","End","Copy")]
  colnames(logReadCountsMode) = c("Chrom","Start","End","log2Ratio")
  write.table(logReadCountsMode,file=paste0(name,".Copywriter.log2RR.Mode.txt"),quote=F,sep="\t",row.names=F,col.names=T)
  
  logReadCounts = logReadCounts[,c("Chromosome","Start","End","Copy")]
  colnames(logReadCounts) = c("Chrom","Start","End","log2Ratio")
  write.table(logReadCounts,file=paste0(name,".Copywriter.log2RR.MAD.txt"),quote=F,sep="\t",row.names=F,col.names=T)
}

