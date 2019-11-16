#module load python/2.7-anaconda-4.4 && source activate Cori_new
# Created by Sofia Medina, Rokhsar lab, UC Berkeley - Nov 15, 2019
# To run: python Exonerate_parser.py -e EXONERATE_Xenla10.1/xl_mgc_cds_nt_0.exonerate -o EXONERATE_Xenla10.1_OUT_TST -a /global/cscratch1/sd/sofiamr/SAPS/X_laevis_V10.1/Genome_Completeness/Fasta_CDS/xl_mgc_cds_nt_anot.tab

import re
import numpy as np
import os
import pandas as pd 
import argparse, os, sys


def main():
    args = parse_args()
    fname = args.exonerate
    annotation_file = args.annotation
    save_output_dir = args.outputDir
    CDS_anot_file   = args.annotation
    
    if save_output_dir.endswith('/')  == False:
        save_output_dir = ''.join((save_output_dir,'/'))
        
    if not os.path.exists(save_output_dir):
        os.makedirs(save_output_dir)
        print "Created directory for:\t\t\n", save_output_dir
        
    out_aln_file_name   =  ''.join((save_output_dir,fname.split('/')[-1].replace('.exonerate','.aln')))
    exonerate_parser(fname, out_aln_file_name, CDS_anot_file, save_output_dir)
    print "Done running program"
    return()
    
    

def parse_args():
    parser = argparse.ArgumentParser(description='Split large fasta file into multiple files of equal number of sequences. The program also create a\
n annotation file that is required for Exonerate (cds2genome)')
    parser.add_argument('-a', '--annotation', metavar='STR', help='CDS Annotation file', type=str, default='')
    parser.add_argument('-o', '--outputDir', metavar='STR', help='Output directory', type=str, default='Exonerate_summary')
    parser.add_argument('-e', '--exonerate', metavar='STR', help='Exonerate output file', type=str, default='')
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    if len(args.exonerate) == 0:
        print "Please provide Exonerate output"
        sys.exit(2)
    if len(args.annotation) == 0:
        print "Please provide CDS Annotation file"
        sys.exit(2)

    return(args)
    
def extract_query_info(line):
    gene_info = line.split('Query: ')[1]#.split(' ')[0]
    if 'lcl' in gene_info:
        gene_id =  gene_info.replace("lcl|","").split("_cds_")[0].replace(' ','')
        prot_id =  gene_info.replace("lcl|","").split("_cds_")[1].split(' ')[0].replace(' ','')
        gene_name = line.split("Query: ")[1].split(' ')[1].split("[gene=")[1].replace(']','').replace(' ','')
    else:
        gene_id = gene_info.split(' ')[0]
        prot_id = ' '.join((gene_info.split(' ')[1:]))
        gene_name =  prot_id
    return(gene_id, prot_id, gene_name)

def identify_genomic_indels(current_line):
    aligment =  ''.join(( str(current_line['Start']),' ', current_line[1],' ',str(current_line['End']),' ',current_line['Strand'],'\n',
str(current_line['Start']),' ', current_line[2],' ',str(current_line['End']),' ',current_line['Strand'],'\n',
str(current_line['Start']),' ', current_line[3],' ',str(current_line['End']),' ',current_line['Strand'],'\n',
str(current_line['Start']),' ', current_line[4],' ',str(current_line['End']),' ',current_line['Strand'],'\n',
str(current_line['Start']),' ', current_line[5],' ',str(current_line['End']),' ',current_line['Strand']))
    
    start        = current_line['Start']
    end          = current_line['End']
    strand       = current_line['Strand']
    tot_aln_len  = len(current_line[5])
    nsplit       = len(list(current_line[5].split("...")))
    if  nsplit >1:
        first_range  = [0,len(current_line[5].split("...")[0])]
        second_range = [tot_aln_len-len(current_line[5].split("...")[-1]),tot_aln_len]
        segments = [current_line[3][first_range[0]:first_range[1]], current_line[3][second_range[0]:second_range[1]]]
        nsegments = 2
    if nsplit==1:
        segments = list([current_line[3]])
        nsegments = 1
    
    genomic_positions_indels = []
    
    for chunk in range(0,nsegments):
        detect_frameshift = segments[chunk]
        #print 'detect_frameshift', detect_frameshift
        pos = np.array([pos for pos, char in enumerate(detect_frameshift) if char == "#"])
        if chunk == 0:
            if strand == '+': genomic_position  =  start + pos;  genomic_positions_indels.extend(genomic_position)
            if strand == '-': genomic_position  =  start - pos;  genomic_positions_indels.extend(genomic_position)
        if chunk == 1:
            distance_from_end =  len(detect_frameshift) - pos
            if strand == '+': genomic_position  =  end - distance_from_end; genomic_positions_indels.extend(genomic_position)
            if strand == '-': genomic_position  =  end + distance_from_end; genomic_positions_indels.extend(genomic_position)
    
    to_return = (sorted(list(set(genomic_positions_indels))),aligment)
    return(to_return)


def exonerate_parser(fname, out_aln_file_name, anot_file, save_output_dir):
    qbfile              = open(fname,"r")
    out_aln_file        = open(out_aln_file_name,'w')
    out_table           = out_aln_file_name.replace('.aln','.tab')
      
        
    dictionary_features = dict()
    first_time          = True
    found_frameshift    = False
    hsp_num             = 0
    ct                  = 0
    flag                = 0
    
    for aline in qbfile.readlines():
        line = aline.replace('\n','')
        if ("Query:" in line) :
            hsp_num=hsp_num+1
            dictionary_features['Indels',hsp_num] = []
            first_time = True
            if "_cds_" in line:
                gene_id, prot_id, gene_name = extract_query_info(line)
            
                dictionary_features['Q_gene_id',hsp_num] = gene_id
                dictionary_features['Q_prot_id',hsp_num] = prot_id
                dictionary_features['Q_gene_name',hsp_num] = gene_name
            else:
                gene_id, prot_id, gene_name = extract_query_info(line)
                dictionary_features['Q_gene_id',hsp_num] = gene_id
                dictionary_features['Q_prot_id',hsp_num] = prot_id
                dictionary_features['Q_gene_name',hsp_num] = gene_name
                
                #dictionary_features['Q_gene_id',hsp_num]=line.split("[")[0]
        if ("Raw score:" in line):
            dictionary_features['Score',hsp_num] = int(line.split(":")[1])
          
        if ("Query range:" in line):
            Q_range = line.replace('Query range:','').replace('->',',').replace(' ','')
            Q_start,Q_end = map(int, Q_range.split(','))
            dictionary_features['Q_start',hsp_num] = Q_start
            dictionary_features['Q_end',hsp_num]   = Q_end
            q_strand = Q_end - Q_start
            if q_strand>0:
                dictionary_features['Q_strand',hsp_num]  = "+"
                q_strand ="+"
            else:
                dictionary_features['Q_strand',hsp_num]  = "-"
                q_strand ='-'
        if ("Target:" in line) :
            target_name = line.split('Target: ')[1].split(' ')[0]
            if 'revcomp' in line:
                dictionary_features['T_strand',hsp_num]  = "-"
            else:
                dictionary_features['T_strand',hsp_num]  = "+"
            dictionary_features['Target',hsp_num] = target_name
        if ("Target range:" in line):
            T_range = line.replace('Target range:','').replace('->',',').replace(' ','')
            T_start,T_end = map(int, T_range.split(','))
            dictionary_features['T_start',hsp_num] = T_start
            dictionary_features['T_end',hsp_num]   = T_end
            T_range = T_range.split(',')
            t_strand = int(T_end)- int(T_start)
            counter_line_aln = -1
    
            flag=0
        if line.startswith('['):
            dictionary_features['Pct_ID',hsp_num] = float(line.split('\t')[0].replace("[","").replace("]",""))
        if 'cigar' in line:  #This means we are done reading hsp
            flag = 1
            
        if flag ==0:
            if line == '': 
                counter_line_aln =0; 
                current_line = dict()
            else: counter_line_aln = counter_line_aln +1
            if (" : " in line):
                if counter_line_aln==5:
                    current_line['Start'] = int(aline.split(':')[0])
                    current_line['End']   = int(aline.split(':')[-1])
                    current_line['Strand'] = dictionary_features['T_strand',hsp_num]
                spaces_before_alingment = line.find(":")+2
                current_line[counter_line_aln] =  line[spaces_before_alingment:-spaces_before_alingment+1]
            else: 
                if counter_line_aln>0:
                    current_line[counter_line_aln] =  line[spaces_before_alingment:]   
            if ("|" in line):
                line_strings = re.sub('\d','', line[spaces_before_alingment:].replace(' bp','').replace('{','').replace('}',''))
                if "#" in line_strings :
                    quept_query_line = line
                    ct=1
                    if "#" in line:
                        found_frameshift = True
        ct = ct +1
        if ct >3: 
            if found_frameshift == True:
                indels, alignment  = identify_genomic_indels(current_line)
                dictionary_features['Indels',hsp_num].extend(indels)
                to_print_aln = ''.join(("Query:",gene_id,",",gene_name,"\t","Indels:\t",target_name, str(indels),'\n',alignment,'\n\n'))
                #print to_print_aln
                out_aln_file.write(to_print_aln)
            found_frameshift=False
    
    out_aln_file.close()
    
    
    ### Making dataframe from dictionary:
    
    annotation = pd.read_csv(anot_file, sep="\t", names=['Q_gene_id','strand','start','length','na']).reset_index()
    annotation = annotation[['Q_gene_id','length']]
    
    Exonerate_summary = pd.DataFrame.from_dict(dictionary_features, orient='index')
    Exonerate_summary['hsp'] = pd.Series(Exonerate_summary.index.tolist()).apply(lambda x: x[1]).tolist()
    Exonerate_summary['Col'] = pd.Series(Exonerate_summary.index.tolist()).apply(lambda x: x[0]).tolist()
    Exonerate_summary = Exonerate_summary[['hsp','Col',0]].sort_values('hsp')
    Exonerate_summary = Exonerate_summary.pivot( columns='Col',index='hsp',values=0)
    Exonerate_summary['Q_len'] = abs(Exonerate_summary.Q_end - Exonerate_summary.Q_start)
    Exonerate_summary['T_len'] = abs(Exonerate_summary.T_end - Exonerate_summary.T_start)
    
    #Merge Dataframe that we just generated with the sequence length to obtain coverage.
    Exonerate_summary = pd.merge(annotation, Exonerate_summary, on='Q_gene_id')
    Exonerate_summary['Q_coverage'] = Exonerate_summary.Q_len/Exonerate_summary.length*100
    Exonerate_summary['Score_rank'] = Exonerate_summary.Score / Exonerate_summary.groupby('Q_gene_id').Score.transform(np.max)
    Exonerate_summary['nIndels'] = Exonerate_summary.Indels.apply(lambda x: len(x))
    cols = ['Q_gene_name','Q_gene_id','length','Q_coverage','Pct_ID','Target','T_start','T_end','T_strand','Score','Score_rank','Indels', 'nIndels']

    print "# genes found: ", Exonerate_summary[Exonerate_summary.Score_rank==1].sort_values(by=['Q_gene_id','Score_rank','Q_coverage'], ascending =False)[cols].shape[0]
    
    Exonerate_summary.sort_values(by=['Q_gene_id','Score_rank','Score','Q_coverage'], ascending =False)[cols].to_csv(out_table, sep="\t", index=False)

    print "Outputs saved as:\n", out_table,'\n', out_aln_file_name
    print "Done with ", fname
    print "Outputs find in:", save_output_dir
    
    Exonerate_summary[Exonerate_summary.Score_rank==1].sort_values(by=['Q_gene_id','Score_rank','Score','Q_coverage'], ascending =False)[cols].head()
    return()
    
    
if __name__=='__main__':
    main()
