#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example Command
# python concat_seqs_order.py --fasta_file_list_infile /Users/kevin.muirhead/Desktop/phytoHybGenes/gene_sequences/hybrid_genes_files.txt --gene_name_list "cpn60,secY,secA,nusA,tuf" --sample_name "BbSP" --organism_name "Candidatus Phytoplasma asteris" --output_dir /Users/kevin.muirhead/Desktop/phytoHybGenes

parser = argparse.ArgumentParser()

fasta_file_list_infile = None
gene_name_list = None
sample_name = None
output_dir = None

parser.add_argument('--fasta_file_list_infile', action='store', dest='fasta_file_list_infile',
                    help='input fasta file path list as input. (i.e. fasta_file_list.txt)')
parser.add_argument('--gene_name_list', action='store', dest='gene_name_list',
                    help='The ordered gene name list (i.e. "cpn60,secY,secA,nusA,tuf")')
parser.add_argument('--sample_name', action='store', dest='sample_name',
                    help='sample_name (i.e. "BbSP")')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_file_list_infile = results.fasta_file_list_infile
gene_name_list = results.gene_name_list
sample_name = results.sample_name
output_dir = results.output_dir

if(fasta_file_list_infile == None):
    print('\n')
    print('error: please use the --fasta_file_list_infile option to specify the input fasta file path list as input')
    print('fasta_file_list_infile =' + ' ' + str(fasta_file_list_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(gene_name_list == None):
    print('\n')
    print('error: please use the --gene_name_list option to the ordered gene name list as input')
    print('gene_name_list =' + ' ' + str(gene_name_list))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(sample_name == None):
    print('\n')
    print('error: please use the --sample_name option to specify the sample_name as input')
    print('sample_name =' + ' ' + str(sample_name))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


fasta_file_list_input_file = open(fasta_file_list_infile, "r")
counter = 0
gene_sequences = {}
for fasta_infile in fasta_file_list_input_file.readlines():
    
    # Need to strip the newline character using the readlines() function for parsing text files.
    fasta_infile = fasta_infile.rstrip('\n')
    
    fasta_basename = os.path.basename(fasta_infile)
    fasta_filename = os.path.splitext(fasta_basename)[0]
   
    print(fasta_filename)
    
    (sample_filename, gene_filename) = fasta_filename.split("__")
    
    for record in SeqIO.parse(fasta_infile, "fasta"):

        seq_id = record.id
        desc = record.description
        seq = record.seq

        #print(seq_id)
        print(gene_filename)
        #print(desc)
        seq_length = len(seq)
        #print(seq_length)
        gene_sequences[gene_filename] = seq

# Concatenate the sequences in order of gene_name i.e "cpn60,secY,secA,nusA,tuf"
concat_seq = ""
concat_length_check = 0
length_check_list = []
genes_present_list =[]
for gene_name in gene_name_list.split(","):
    #print(gene_name)
    if(gene_name in gene_sequences):
        concat_seq += str(gene_sequences[gene_name])
    
        seq_length = len(str(gene_sequences[gene_name]))
        #print(str(seq_length))
        genes_present_list.append(gene_name)
    else:

        concat_seq += ""
        seq_length = len("")

    length_check_list.append(str(seq_length))
    
    concat_length_check += seq_length

# Print the concatenated hybrid gene fasta file.
fasta_outfile = os.path.join(output_dir, "_".join([sample_name,"concat_hybrid_genes.fasta"]))
fasta_output_file = open(fasta_outfile, "w+")
#print(concat_seq)

# Get the length of the concatenated sequence.
concat_seq_length = len(str(concat_seq))

print("This is is a length check for the calculation of the concatenated length.")
print(gene_name_list)
print(" + ".join(length_check_list) + " = " + str(concat_seq_length))
print(concat_seq_length)

# The header description
desc = " ".join(["\"" + ",".join(genes_present_list) + "\"","length=" + str(concat_seq_length) + "bp"])
#print(desc)

# Make a new sequence record object in fasta format with header and concatenated sequence.
concat_record = SeqRecord(
    Seq(concat_seq),
    id=sample_name,
    name="",
    description=desc
)

SeqIO.write(concat_record, fasta_output_file, "fasta")
fasta_output_file.close()


