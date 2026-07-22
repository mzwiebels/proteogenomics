#!/usr/bin/Rscript

##########################################################################################
##
## SNV_SelectOutputSS.R
##
## Filters output for annotated files (single-sample).
##
##########################################################################################

message("\n##Filter SNV Output###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-t", "--type"),type="character",default=NULL,help="sample type (Tumor | Normal)"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-c", "--CGC"),type="character",default=NULL,help="CGC annotation file"),
  make_option(c("-r", "--TruSight"),type="character",default=NULL,help="TruSight annotation file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$name)){
  print_help(opt_parser)
  stop("You have to specify a sample name.",call.=FALSE)
} 
if(is.null(opt$type) | !opt$type %in% c("Tumor","Normal")){
  print_help(opt_parser)
  stop("You have to specify the input type (Tumor | Normal).",call.=FALSE)
}
if(is.null(opt$species) | !opt$species %in% c("Human","Mouse")){
  print_help(opt_parser)
  stop("You have to specify the species (Human | Mouse).",call.=FALSE)
}


# main program
if(!require("data.table")) install.packages("data.table")
suppressMessages(library(data.table))

name = opt$name
type = opt$type
species= opt$species
CGC = opt$CGC
TruSight = opt$TruSight


file=as.data.frame(fread(paste(name,"/results/Mutect2/",name,".",type,".Mutect2.txt",sep="")),header=T,sep="\t")

if (species=="Human"){
	file$AF[is.na(file$AF)] = 0
	file$AC[is.na(file$AC)] = 0
	file$AN[is.na(file$AN)] = 0
  
	sel = file[file[,"AF"] <0.1 & file[,"G5"]=="FALSE" & file[,"AN"] < 100 | file[,"AF"] <0.01 & file[,"AN"] >= 100 & file[,"G5"]=="FALSE",] 
	sel[is.na(sel)] = " "
	
	write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
  
	
	sel = sel[sel[,"ANN[*].IMPACT"] %in% c("HIGH", "MODERATE"),]
	sel[is.na(sel)] = " "
	
	write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.OnlyImpact.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
  
	if(!(is.null(CGC) | CGC=="NULL")){
	  CGC = read.delim(CGC,header=T,sep="\t")
	  sel = sel[sel[,"ANN[*].GENE"]%in% CGC[,1],]
	  sel[is.na(sel)] = " "
	  
	  write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
	}
	
	if(!(is.null(TruSight) | TruSight=="NULL")){
  	TruSight = read.delim(TruSight,header=T,sep="\t")
  	sel = sel[sel[,"ANN[*].GENE"]%in% TruSight[,1],]
  	sel[is.na(sel)] = " "
  	
  	write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
	}
} else if (species=="Mouse"){
	sel = file[file[,"ANN[*].IMPACT"] %in% c("HIGH", "MODERATE"),]
	sel[is.na(sel)] = " "
	
	write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.OnlyImpact.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
  
	if(!(is.null(CGC) | CGC=="NULL")){
  	CGC = read.delim(CGC,header=T,sep="\t",fill=T)
  	sel = sel[sel[,"ANN[*].GENE"]%in% CGC[,1],]
  	sel[is.na(sel)] = " "
  	
  	write.table(sel,paste(name,"/results/Mutect2/",name,".",type,".Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
	}
}