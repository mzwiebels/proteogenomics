#!/usr/bin/Rscript

##########################################################################################
##
## CNV_RunCopywriter.R
##
## Run Copywriter on matched tumour-normal .bam-files using 20kb windows.
##
##########################################################################################
message("\n##CopywriteR CNA estimation###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-t", "--threads"),type="integer",default=8,help="number of threads [default: 8]"),
  make_option(c("-m", "--runmode"),type="character",default=NULL,help="runmode (MS | SS)"),
  make_option(c("-g", "--genome_dir"),type="character",default=NULL,help="path to reference data directory"),
  make_option(c("-c", "--centromere_file"),type="character",default=NULL,help="path to centromere exclusion file"),
  make_option(c("-v", "--varregion_file"),type="character",default=NULL,help="path to reference genome fasta file"),
  make_option(c("-r", "--resolution"),type="integer",default=20000,help="genomic resolution to call copy number alterations [default: 20,000]"),
  make_option(c("-p", "--type"),type="character",default=NULL,help="sample type (Tumor | Normal)")
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
if(is.null(opt$runmode) | !opt$runmode %in% c("MS","SS")){
  print_help(opt_parser)
  stop("You have to specify the species (Human | Mouse).",call.=FALSE)
} 


# main program
if(!require("CopywriteR")) BiocManager::install("CopywriteR")
if(!require("GenomeInfoDb")) BiocManager::install("GenomeInfoDb")
if(!require("GenomicRanges")) install.packages("GenomicRanges")
if(!require("naturalsort")) BiocManager::install("naturalsort")
suppressMessages(library(CopywriteR))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(GenomicRanges))
suppressMessages(library(naturalsort))

name = opt$name
species = opt$species
threads = opt$threads
runmode = opt$runmode
genome_dir = opt$genome_dir
centromere_file = opt$centromere_file
varregions_file = opt$varregions_file
resolution = opt$resolution
type = opt$type


tumor_bam = paste0(name,"/results/bam/",name,".Tumor.bam")
normal_bam = paste0(name,"/results/bam/",name,".Normal.bam")

if (runmode == "MS"){
	sample.control = data.frame(samples=c(normal_bam,tumor_bam),controls=c(normal_bam,normal_bam))
} else if (runmode == "SS"){
	if(type == "Tumor"){
	  sample.control = data.frame(samples=c(tumor_bam),controls=c(tumor_bam))
	} else if(type == "Normal"){
	  sample.control = data.frame(samples=c(normal_bam),controls=c(normal_bam))
	}
}

resolution = resolution/1000
if(species == "Human"){
  reference_files = paste0(genome_dir,"/hg38_",resolution,"kb")
} else if(species == "Mouse"){
  reference_files = paste0(genome_dir,"/mm10_",resolution,"kb")
}

bp.param = SnowParam(workers=threads,type="SOCK")

CopywriteR(sample.control = sample.control,
             destination.folder = file.path(paste0(name,"/results/Copywriter")),
             reference.folder = file.path(reference_files),
             bp.param = bp.param)

log2.reads = read.table(paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"),header=T,sep="\t",check.names=FALSE)
file.copy(paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"),paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts_backup.igv"),overwrite=T)
if (runmode == "MS"){
  log2.reads.GR = GRanges(log2.reads$Chromosome,IRanges(log2.reads$Start,log2.reads$End),Feature=as.character(log2.reads$Feature),Normal=log2.reads[,5],Tumor=log2.reads[,6])
} else if (runmode == "SS"){
  if(type == "Tumor"){
    log2.reads.GR = GRanges(log2.reads$Chromosome,IRanges(log2.reads$Start,log2.reads$End),Feature=as.character(log2.reads$Feature),Tumor=log2.reads[,5])
  } else if(type == "Normal"){
    log2.reads.GR = GRanges(log2.reads$Chromosome,IRanges(log2.reads$Start,log2.reads$End),Feature=as.character(log2.reads$Feature),Normal=log2.reads[,5])
  }
}


# remove regions with increased variability for mice and centromere regions for humans
if(species == "Human"){
	filter = read.table(centromere_file,header=F,sep="\t")
	flankLength = 5000000
}
if(species == "Mouse"){
  filter = read.table(varregions_file,header=F,sep="\t")
	flankLength = 0
}

colnames(filter)[1:3] <- c("Chromosome","Start","End")
filter$Start <- filter$Start - flankLength
filter$End <- filter$End + flankLength
filter = GRanges(filter$Chromosome,IRanges(filter$Start,filter$End))

hits <- findOverlaps(query=log2.reads.GR,subject=filter)
ind <- queryHits(hits)
log2.reads.GR=(log2.reads.GR[-ind, ])
message("Removed ", length(ind), " bins near centromeres.")

log2.reads.fixed = as.data.frame(log2.reads.GR)
if (runmode == "MS"){
  log2.reads.fixed = log2.reads.fixed[,c("seqnames","start","end","Feature","Normal","Tumor")]
  colnames(log2.reads.fixed) = c("Chromosome", "Start", "End", "Feature",paste0("log2.",name,".Normal.bam"),paste0("log2.",name,".Tumor.bam"))
} else if (runmode == "SS"){
  if(type == "Tumor"){
    log2.reads.fixed = log2.reads.fixed[,c("seqnames","start","end","Feature","Tumor")]
    colnames(log2.reads.fixed) = c("Chromosome", "Start", "End", "Feature",paste0("log2.",name,".Tumor.bam"))
  } else if(type == "Normal"){
    log2.reads.fixed = log2.reads.fixed[,c("seqnames","start","end","Feature","Normal")]
    colnames(log2.reads.fixed) = c("Chromosome", "Start", "End", "Feature",paste0("log2.",name,".Normal.bam"))
  }
}
write.table(log2.reads.fixed,paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"),sep="\t",quote=F,row.names=F,col.names=T)
plotCNA(destination.folder = file.path(paste0(name,"/results/Copywriter")),min.width=5,alpha=0.001,undo.splits="sdundo",undo.SD=2)