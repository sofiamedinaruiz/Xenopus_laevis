#module load python/2.7-anaconda-4.4 && source activate Cori_new
# Created by Sofia Medina, Rokhsar lab, UC Berkeley - Nov 15, 2019
# Prepare the input files for Exonerate (cds2genome). 

import argparse, os, pysam, sys
import pandas as pd
import numpy as np
import os
from Bio import SeqIO


def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print "Created directory:", directory
    return()

def parse_args():
    parser = argparse.ArgumentParser(description='Split large fasta file into multiple files of equal number of sequences. The program also create an annotation file that is required for Exonerate (cds2genome)')
    parser.add_argument('-s', '--split', metavar='INT', help='number of sequences per file', type=int, default=100)
    parser.add_argument('-o', '--output', metavar='STR', help='Name of output directory', type=str, default='Fasta_out')
    parser.add_argument('-d', '--data', metavar='STR', help='Data file from MGCs', type=str, default='')
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
    required = parser.add_argument_group('required arguments')
    required.add_argument("fasta", help="Fasta file to obtain annotation (ideally complete CDS sequences)", type=str)
    args = parser.parse_args()
    
    if (args.fasta.endswith('.fa') | args.fasta.endswith('.fasta')):
        print "Input fasta:", args.fasta
    else:
        print "Please provide fasta file"
        sys.exit(2) 
    return(args)


def main():
    args = parse_args()
    
    out_dir = ''.join((args.output,'/')).replace('//','/')
    ensure_dir(out_dir)
    prefix_fasta = '.'.join((args.fasta.split('.')[:-1]))
    prefix_all = ''.join((out_dir,prefix_fasta,'_'))
    annotation_file = ''.join((out_dir,prefix_fasta,'_anot.tab'))
    
    
    if len(args.data)>0:
        print "Reading", args.data
        mRNA_DATA = pd.read_csv(args.data, sep='\t')
        genbank_dict = mRNA_DATA[['GenBank accession','GenBank def line']].set_index('GenBank accession').to_dict(orient='index')
    
    counter =0
    
    fasta_out_name = ''.join((prefix_all,str(counter),'.fa'))
    
    if os.path.exists(fasta_out_name):
        print "WARNING: FILES ALREADY CREATED: Please erase files or change output directory"
        sys.exit(2) 
        
    fasta_out = open(fasta_out_name, 'w')
    anot_out = open(annotation_file, 'w')
    print "Will save: ", annotation_file
    print "Will save: ", ''.join((prefix_all,str(counter),'.fa'))

    for record in SeqIO.parse(args.fasta, "fasta"):
        if counter < 100000000:
            anot = '\t'.join((record.id,'+',str(1),str(len(record.seq))))
            anot_out.write((''.join((anot,'\n')))) 
            sequence  = ''.join(('>',record.id,' ',str(genbank_dict[record.id]['GenBank def line']),'\n',str(record.seq),'\n'))
            fasta_out.write(sequence)
            if np.mod(counter, args.split) == 0:
                fasta_out.close()
                fasta_out_name = ''.join((prefix_all,str(counter),'.fa'))
                print "Will save: ", fasta_out_name
                fasta_out = open(''.join((prefix_all,str(counter),'.fa')), 'w')

            counter=counter+1
    
    anot_out.close()
    fasta_out.close()
    
    
if __name__=='__main__':
    main()
