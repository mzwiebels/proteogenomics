#!/usr/bin/Rscript

##########################################################################################
##
## LOH_GenerateVariantTable.R
##
## Calculate datapoints needed for plotting LOH while transforming them to B-allele frequencies.
##
##########################################################################################
message("\n##Filter SNV Output###")
options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-s", "--species"),type="character",default=NULL,help="sample species (Human | Mouse)"),
  make_option(c("-f", "--fasta"),type="character",default=NULL,help="path to genome fasta file"),
  make_option(c("-r", "--repository"),type="character",default="repository",help="path to script repository [default: ./repository/]")
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
suppressMessages(library(data.table))

name = opt$name # used for naming in- and output files
species = opt$species # sample species, 'Mouse' or 'Human'
genome_fasta = opt$fasta # location of the reference genome fasta
repository_dir = opt$repository # location of repository

source(paste0(repository_dir,"/LOH_Library.R"))


#read input files
tumor = as.data.frame(fread(paste0(name,"/results/Mutect2/",name,".Tumor.Mutect2.Positions.txt"), header=T, sep="\t"))
normal = as.data.frame(fread(paste0(name,"/results/Mutect2/",name,".Normal.Mutect2.Positions.txt"), header=T, sep="\t"))

#adjust column names
colnames(tumor) = c("Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount", "MapQ", "BaseQ")
colnames(normal) = c("Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount", "MapQ", "BaseQ")

#filter reads
tumor = LOH_FilterReads(tumor,species)
normal = LOH_FilterReads(normal,species)

#merge both files using UniquePos, reduce the input to heterozygous positions
variants = LOH_MergeVariants(tumor,normal)

#write out table used for plotting of germline variants
write.table(variants,file=paste0(name,"/results/LOH/",name,".VariantsForLOHGermline.txt"),quote=F,sep="\t",row.names=F,col.names=T)

#filter for informative variants (heterozygous in germline)
variants = variants[variants[,"Normal_Freq"] <= 0.7 & variants[,"Normal_Freq"] >= 0.3,]

#define the dictionary, which defines whether an allele im AMBigous or UNAMBigous.
dicts = LOH_DefineDictionaries()

dict_for = dicts[["dict_for"]]
dict_rev = dicts[["dict_rev"]]
dict_unamb = dicts[["dict_unamb"]]

#determine which variants can be assigned A/B immediately (unamb), and which need further lookup (amb)
results = LOH_AssignStatus(variants,dict_unamb)

amb = results[["amb"]]
indel = results[["indel"]]
unamb = results[["unamb"]]

#determine which allele serves as A- and B-allele
LOH_GenerateLookupTable(name,amb)

#extract positions surrounding the variants from the reference genome. These are needed for the assignemnt of the A- and B-Allele.
LOH_ExtractSurroundingNucleotides(name,genome_fasta)

#Re-Import all nucleotides, which surround the requested variants
import = LOH_ImportSurroundingNucleotides(name)

#Using the information from the surrounding nucleotides, all variants in amb are assigned to A or B
amb = LOH_AssignAlleles(amb,import,dict_unamb)

#merge the final tables and export it for plotting
final_variants = rbind(amb,unamb,indel)
write.table(final_variants,file=paste0(name,"/results/LOH/",name,".VariantsForLOH.txt"),quote=F,sep="\t",row.names=F,col.names=T)