#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

fasta_infile = None
min_seq_length = None
fasta_outfile = None

parser.add_argument('-i', action='store', dest='fasta_infile',
                    help='input fasta file path as input. (i.e. sequences.fasta)')
parser.add_argument('-l', action='store', dest='min_seq_length',
                    help='input minimum sequence length. Default: 500')
parser.add_argument('-o', action='store', dest='fasta_outfile',
                    help='output fasta file path as input. (i.e. new_sequences.fasta)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
min_seq_length = results.min_seq_length
fasta_outfile = results.fasta_outfile

if(fasta_infile == None):
	print('\n')
	print('error: please use the -i option to specify the input fasta file path as input')
	print('fasta_infile =' + ' ' + str(fasta_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(min_seq_length == None):
    print('\n')
    print('error: please use the -l option to specify the input minimum sequence length. Default: 500')
    print('fasta_infile =' + ' ' + str(min_seq_length))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(fasta_outfile == None):
    print('\n')
    print('error: please use the -o option to specify the output fasta file path as input')
    print('fasta_outfile =' + ' ' + str(fasta_outfile))
    print('\n')
    parser.print_help()
    sys.exit(1)


basename = os.path.basename(fasta_outfile)
filename = os.path.splitext(basename)[0]
output_dir = os.path.dirname(fasta_outfile)

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

fasta_output_file = open(fasta_outfile, "w+")
for record in SeqIO.parse(fasta_infile, "fasta"):
    seq_id = record.id
    desc = record.description
    seq = record.seq
    seq_length = len(seq)
    
    if(int(seq_length) >= int(min_seq_length)):
    
        SeqIO.write(record, fasta_output_file, "fasta")
        
fasta_output_file.close()


