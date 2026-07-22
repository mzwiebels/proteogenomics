#!/usr/bin/env python
######################################################################################
# 
# VCF Processing Script Version 3 - better handling of the NCBI database parsing
#
######################################################################################

import argparse
import pandas as pd
import re
import subprocess
#import sys
from Bio import SeqIO, Entrez
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

Entrez.email = 'zwiebel@biochem.mpg.de'   # Entrez, the official API for databank queries on NCBI, requires the user to input their email address 

parser = argparse.ArgumentParser(description="vcf_to_fasta converts mutation calls from NGS in vcf format into mutated protein sequences in fasta format for use in proteomics.")
parser.add_argument('-i', '--inputFile', type=str, required=True, help='path to input vcf file')
parser.add_argument('-o', '--outputFile', type=str, default='mutated_peptide_sequences.fasta', help='path to output fasta file')

args = parser.parse_args()

### Function to introduce genetic mutations in mRNA
def create_mutation(gene, mutation):
    gene = MutableSeq(gene)

    if '>' in mutation: # Takes care of all simple substitutions ala C > T, A > G or similar
        match = re.search(r'c.(\d+)([A-Z])>([A-Z])', mutation)
        pos = int(match.group(1))
        ref = match.group(2)
        alt = match.group(3)

        if gene[pos-1] == str(ref):
            gene[pos-1] = alt

        return gene, pos
    else: 
        del_seq = ins_seq = dup_seq = inv_seq = None
        
        if '_' in mutation: # when the mutation spans multiple positions, extract start and stop
            match = re.search(r'c.(\d+)_(\d+)', mutation)
            start, stop = int(match.group(1)), int(match.group(2))
        else:
            match = re.search(r'c.(\d+)', mutation) # otherwise, extract the singular position from the mutation handle and set start and stop to it
            start = stop = int(match.group(1))
                
        new_sequence = gene[:start]
        
        if 'del' in mutation:
            match = re.search(r'del([A-Z]+)', mutation)
            if match:
                del_seq = match.group(1)
            else:
                del_seq = ''
                            
        if 'ins' in mutation:
            match = re.search(r'ins([A-Z]+)', mutation)
            if match:
                ins_seq = match.group(1)
            else:
                ins_seq = ''
                
        if 'dup' in mutation:
            match = re.search(r'dup([A-Z]+)', mutation)
            if match:
                dup_seq = match.group(1)
            else:
                dup_seq = gene[start-1:stop]
                
        if 'inv' in mutation:
            match = re.search(r'inv([A-Z]+)', mutation)
            if match:
                inv_seq = match.group(1)
            else:
                inv_seq = ''
                
        if del_seq is not None and ins_seq is not None:
            new_sequence = new_sequence[:start-1] + ins_seq + gene[stop:]
        elif del_seq is not None:
            new_sequence = new_sequence[:start-1] + gene[stop:]
        elif ins_seq is not None:
            new_sequence = new_sequence + ins_seq + gene[start:]
        elif dup_seq is not None:
            new_sequence = new_sequence[:start-1] + dup_seq + gene[start-1:]
        elif inv_seq is not None:
            new_sequence = new_sequence[:start-1] + inv_seq + gene[start-1:]

        return new_sequence, start


### Stage 1: Loading of the VCF file containing the relavent information 

## Read VCF file into dataframe, sort by gene name
vcf = pd.read_csv(args.inputFile, keep_default_na=False, sep='\t')
vcf = vcf.sort_values(by=['Hugo_Symbol'])

## Filter out all variant types which cannot be processed
patterns_to_filter = ['start_lost','5_prime_UTR_premature_start_codon_gain_variant','intron_variant','intergenic_region','intragenic_variant','stop_retained','synonymous_variant','upstream_gene_variant', 'downstream_gene_variant','3_prime_UTR_variant','5_prime_UTR_variant','splice_donor_variant','splice_acceptor_variant','splice_region_variant','non_coding_transcript_exon_variant']
pattern = '|'.join(patterns_to_filter)  # Create a regex pattern
vcf = vcf[~vcf['Consequence'].str.contains(pattern, regex=True)]
vcf = vcf[~vcf['HGVSc'].astype(str).str.contains('+', regex=False)]

## Reformat short description for protein changes
vcf['HGVSp_Short'] = vcf['HGVSp_Short'].str[2:].str.rstrip('*?')

## Group data frame into a format which contains mutations for each gene as lists wihtin the cells of the data frame
vcf = vcf.groupby(['Hugo_Symbol', 'RefSeq', 'SWISSPROT', 'TREMBL']).agg({'HGVSc': lambda x: list(x), 'HGVSp': lambda x: list(x), 'HGVSp_Short': lambda x: list(x), 'Consequence':lambda x: list(x)}).reset_index()
vcf = vcf.rename(columns={'Hugo_Symbol':'Gene Name', 'RefSeq':'Gene ID', 'SWISSPROT':'Swissprot ID', 'TREMBL':'Trembl ID', 'HGVSc':'Mut Gene', 'HGVSp':'Mut Prot', 'HGVSp_Short':'Mut Prot Short', 'Consequence':'Mut Type'})
vcf['Gene ID'] = vcf['Gene ID'].str.split(',').str[0]
vcf = vcf[vcf['Gene ID'] != '']

## Print little overview of file
list_of_genes = vcf['Gene ID'].to_list()
print(f'\n{len(list_of_genes)} relevant entries have been found within the VCF file and will be processed.\n')

### Stage 2: Downloading entries to the relavent genes for later processing. The database entries are downloaded in batches of 100 and processed later.

## Split up the list of genes into packages of 100
list_of_genes = vcf['Gene ID'].to_list()
batch_size = 100

sequence_records = []

for start in range(0, len(list_of_genes), batch_size):
    end = min(start + batch_size, len(list_of_genes))
    batch_ids = list_of_genes[start:end]

    print(f'Fetching records {start+1} to {end}')
    
    try:
        stream = Entrez.efetch(db='nucleotide', id=','.join(batch_ids), rettype='gb', retmode='text')
        record = list(SeqIO.parse(stream, 'gb'))
        stream.close()

    except Exception as e:
        print(f'The following error occured: {e}')

    sequence_records.extend(record)

print('')
#for record in sequence_records:
#    print(record.description)


### Stage 3: Parsing the downloaded information and creating a FASTA file containing the original proteins and the mutated versions

fasta_records = []

for record in sequence_records:

    # Search for the record ID in the vcf data frame
    current = vcf[vcf['Gene ID'] == record.id]
    
    # Extract the coding and complete sequence as well as gene name from the record list
    sequence = record.seq
    
    for feature in record.features:
        if feature.type == 'CDS':
            coding_start = feature.location.start
        elif feature.type == "gene":
            gene_seq = feature.extract(record).seq
    
    gene_seq = gene_seq[coding_start:]
    gene_seq = gene_seq[:len(gene_seq)//3*3]
    
    # Extract the mutation strings and type for the creation of a mutated sequence
    gene_name = current['Gene Name'].values[0]
    protein_id = current['Swissprot ID'].values[0]
    if protein_id == '':
        protein_id = current['Trembl ID'].values[0]
    mutation_Gene = current['Mut Gene'].values[0]
    mutation_Prot = current['Mut Prot'].values[0]
    mutation_Prot_Short = current['Mut Prot Short'].values[0]
    mutation_Type = current['Mut Type'].values[0]

    # First create an entry into the FASTA file for the original protein
    #description_original = 'Original of ' + record.name + ', ' + record.description
    #current_seqrecord = SeqRecord(Seq(gene_seq.translate(to_stop=True)), id=record.name, description=description_original)
    #fasta_records.append(current_seqrecord)

    # Create a muation for the sequence for each mutation string found in mutG
    for mutG, mutP, mutPS, mutT in zip(mutation_Gene, mutation_Prot, mutation_Prot_Short, mutation_Type):
        try:
            mutated_sequence, pos_of_mut = create_mutation(gene_seq, mutG)
        except Exception as e:
            print(f'Exception found on gene {gene_name} mutation {mutG}')
            print(e)

        try:
            if int(subprocess.run(['grep', '-c', protein_id, '/fs/pool/pool-mann-projects/MaxZ/MoCaSeq/data/fasta_headers.txt'], capture_output=True, text=True).stdout.strip()) == 1:
                fasta_header = subprocess.run(['grep', protein_id, '/fs/pool/pool-mann-projects/MaxZ/MoCaSeq/data/fasta_headers.txt'], capture_output=True, text=True).stdout.strip()
            else:
                if int(subprocess.run(['grep', '-c', f'sp.*GN={gene_name} ', '/fs/pool/pool-mann-projects/MaxZ/MoCaSeq/data/fasta_headers.txt'], capture_output=True, text=True).stdout.strip()) == 1:
                    fasta_header = subprocess.run(['grep', f'sp.*GN={gene_name} ', '/fs/pool/pool-mann-projects/MaxZ/MoCaSeq/data/fasta_headers.txt'], capture_output=True, text=True).stdout.strip()
                else:
                    fasta_header = subprocess.run(['grep', '-m1', f'tr.*GN={gene_name} ', '/fs/pool/pool-mann-projects/MaxZ/MoCaSeq/data/fasta_headers.txt'], capture_output=True, text=True).stdout.strip()
            
            # Parse the fasta header correctly
            parts = fasta_header.split(' ', 1)
            name_parts = parts[0].split('|', 2)
            description_parts = parts[1].split('OS=', 1)
            
            custom_id = f"{name_parts[0]}|{name_parts[1]}_{mutPS}|{name_parts[2]}"
            custom_description = f"{description_parts[0]}{mutP} OS={description_parts[1]}"
            
            mutation_seqrecord = SeqRecord(Seq(mutated_sequence.translate(to_stop=True)), id=custom_id, description=custom_description)
            fasta_records.append(mutation_seqrecord)
        except Exception as e:
            print(f'{gene_name} could not be associated to an existing Swissprot entry. Error: {str(e)}')

        # Test-Block 
        #ind = pos_of_mut//3
        #with open('output.txt', 'a') as f:
        #    sys.stdout = f
        #    print(len(gene_seq))
        #    print(f'Mutated sequence for {mutG} in {record.name}')
        #    print(f'Original sequence starting at {pos_of_mut-6}')
        #    print(f'{gene_seq[pos_of_mut-7:pos_of_mut-1]} {gene_seq[pos_of_mut-1:pos_of_mut]} {gene_seq[pos_of_mut:pos_of_mut+6]}')
        #    print(f'Mutated sequence starting at {pos_of_mut-6}')
        #    print(f'{mutated_sequence[pos_of_mut-7:pos_of_mut-1]} {mutated_sequence[pos_of_mut-1:pos_of_mut]} {mutated_sequence[pos_of_mut:pos_of_mut+6]}')
        #    print(f'Sequence on protein level with {mutP}, {mutT}')
        #    print(f'Original translation starting at {ind-6}')
        #    print(f'{gene_seq.translate()[ind-7:ind-1]} {gene_seq.translate()[ind-1:ind+30]}')
        #    print(f'Mutated translation starting at {ind-6}')
        #    print(f'{mutated_sequence.translate()[ind-7:ind-1]} {mutated_sequence.translate()[ind-1:ind+30]}')
        #    print('-'*20)
        #sys.stdout = sys.__stdout__  # Restore stdout
        

SeqIO.write(fasta_records, args.outputFile, 'fasta')
