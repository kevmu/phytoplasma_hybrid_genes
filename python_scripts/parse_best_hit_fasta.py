#!/usr/bin/python
import os
import sys
import re
import csv
import natsort
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


### Example Command
#pip install natsort

parser = argparse.ArgumentParser()

fasta_infile = None
blast_results_infile = None
output_dir = None

parser.add_argument('--fasta_infile', action='store', dest='fasta_infile',
                    help='input fasta file as input. (i.e. fasta_file.fasta)')
parser.add_argument('--blast_results_infile', action='store', dest='blast_results_infile',
                    help='blast results input file as input (i.e. *_blastn.tsv)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
blast_results_infile = results.blast_results_infile
output_dir = results.output_dir

if(fasta_infile == None):
    print('\n')
    print('error: please use the --fasta_infile option to specify the input fasta file as input')
    print('fasta_infile =' + ' ' + str(fasta_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(blast_results_infile == None):
    print('\n')
    print('error: please use the --blast_results_infile option to specify the blast results input file as input')
    print('blast_results_infile =' + ' ' + str(blast_results_infile))
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

blast_results = {}
i = 0
with open(blast_results_infile, "r") as blast_results_input_file:
    csv_reader = csv.reader(blast_results_input_file, delimiter='\t')
    for row in csv_reader:
        if(i != 0):
            
            (qseqid,qacc,qlen,sseqid,sacc,slen,stitle,qstart,qend,sstart,send,length,evalue,bitscore,pident,qcovs,qcovhsp,nident,positive,mismatch,gaps,qframe,sframe) = row
            
            if(not(bitscore in blast_results)):
                blast_results[bitscore] = []
            if(bitscore in blast_results):
                blast_results[bitscore].append([qseqid,qacc,qlen,sseqid,sacc,slen,stitle,qstart,qend,sstart,send,length,evalue,bitscore,pident,qcovs,qcovhsp,nident,positive,mismatch,gaps,qframe,sframe])
                    
        i = i + 1

fasta_record_dict = {}
for record in SeqIO.parse(fasta_infile, "fasta"):

    seq_id = record.id
    desc = record.description
    seq = record.seq

    seq_length = len(seq)
    print(seq_id)
    print(seq)
    fasta_record_dict[seq_id] = record

basename = os.path.basename(fasta_infile)
filename = os.path.splitext(basename)[0]

fasta_outfile = os.path.join(output_dir, filename + "_best_hit.fasta")
fasta_output_file = open(fasta_outfile, "w+")

# Grab the best hit within the blast results based on bitscore.
bitscore = natsort.natsorted(set(blast_results.keys()))[-1]
for row in blast_results[bitscore]:
    print(row)
    (qseqid,qacc,qlen,sseqid,sacc,slen,stitle,qstart,qend,sstart,send,length,evalue,bitscore,pident,qcovs,qcovhsp,nident,positive,mismatch,gaps,qframe,sframe) = row
    
    record = fasta_record_dict[qseqid]
    
    if(int(qframe) == 1):
        SeqIO.write(record, fasta_output_file, "fasta")
    if(int(qframe) == -1):
        record.seq = str(record.seq.reverse_complement())
        SeqIO.write(record, fasta_output_file, "fasta")
        
fasta_output_file.close()
