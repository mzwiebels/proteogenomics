#!/usr/bin/Rscript

##########################################################################################
##
## Preparation_AnnotateGencode.R
##
## Creates .rds files for Gencode gene and exome annotations
##
##########################################################################################
args <- commandArgs(TRUE)

infile = args[1]
version = args[2]
gencodeVersion = args[3]

gencode = read.csv(infile,header=F,sep="\t")
gencode = gencode[,c(1,3,4,5,7,9)]
colnames(gencode) = c("chr","class","start","end","strand","annotation")

genes = gencode[gencode$class=="gene",]
gene_annotations = t(sapply(genes$annotation, function(x){
  tmp = unlist(strsplit(x,"; "))[1:3]
  sapply(strsplit(tmp," "), function(y) return(y[2]))
}))
colnames(gene_annotations) = c("geneIDVersion","geneType","geneName")
genes = cbind(genes,gene_annotations)
genes = genes[,-6]
genes$geneID = sapply(genes$geneIDVersion, function(x) unlist(strsplit(x,"\\."))[1])
saveRDS(genes,file=paste0("ref/",version,"/",version,".",gencodeVersion,".Genes.rds"))

exons = gencode[gencode$class=="exon",]
exon_annotations = t(sapply(exons$annotation, function(x){
  tmp = unlist(strsplit(x,"; "))[c(1:5,7)]
  sapply(strsplit(tmp," "), function(y) return(y[2]))
}))
colnames(exon_annotations) = c("geneID","transcriptID","geneType","geneName","transcriptType","exonNum")
exons = cbind(exons,exon_annotations)
exons = exons[,-6]
saveRDS(exons,file=paste0("ref/",version,"/",version,".",gencodeVersion,".Exons.rds"))
