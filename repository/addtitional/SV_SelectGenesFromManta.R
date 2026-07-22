#!/usr/bin/Rscript

##########################################################################################
##
## SV_SelectGenesFromManta.R
##
## Selects genes from Manta for use in oncoprints.
##
##########################################################################################

message("\n##Breakpoint Cluster Analysis###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-t", "--type"),type="character",default="",help="sample type (Tumor | Normal), or none. Defaults to none ('')")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$name)){
  print_help(opt_parser)
  stop("You have to the sample name.",call.=FALSE)
}
if(is.null(opt$type)){
  print_help(opt_parser)
  stop("You have to specify the sample type (Tumor | Normal | '').",call.=FALSE)
}

# main program
if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")
f(!require("tidyr")) install.packages("dplyr")
f(!require("readr")) install.packages("dplyr")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))


if (type==""){
	tab = read_tsv(paste(name,"/results/Manta/",name,".Manta.txt",sep=""))
} else{
	tab = read_tsv(paste(name,"/results/Manta/",name,".",type,".Manta.txt",sep=""))
}

colnames(tab) = gsub("\\[\\*]","",colnames(tab))
colnames(tab) = gsub("\\[","",colnames(tab))
colnames(tab) = gsub("\\]","",colnames(tab))

tab = tab %>% 
  mutate_if(is.integer, as.numeric) %>%
  mutate(GENTumor.SR.AF=GENTumor.SR1/(GENTumor.SR1+GENTumor.SR0)) %>%
  mutate(GENTumor.PR.AF=GENTumor.PR1/(GENTumor.PR1+GENTumor.PR0)) %>%
  filter(GENTumor.SR1 + GENTumor.SR0 >= 10 | GENTumor.PR1 + GENTumor.PR0 >= 10) %>%
  filter(ANN.IMPACT %in% c("HIGH","MODERATE")) %>%
  filter(grepl("fusion|frameshift|stop|splice|start",ANN.EFFECT)) %>%
  separate(ANN.GENE,into=c("Gene1","Gene2"),sep="&") %>% 
  select(c("Gene1","Gene2","ANN.EFFECT")) %>%
  gather(key,Gene,1:2) %>%
  select(-key) %>%
  filter(!is.na(Gene)) %>% 
  separate(ANN.EFFECT,into=c("EFF1","EFF2"),sep="&") %>%
  gather(key,Effect,1:2) %>% 
  select(-key) %>%
  filter(!is.na(Effect)) %>%
  distinct()

if (type==""){
	write.table(tab,paste(name,"/results/Manta/",name,".Manta.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
} else{
	write.table(tab,paste(name,"/results/Manta/",name,".",type,".Manta.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}
