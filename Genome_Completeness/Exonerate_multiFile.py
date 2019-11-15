#module load python/2.7-anaconda-4.4 && source activate Cori_new
# Created by Sofia Medina, Rokhsar lab, UC Berkeley - Nov 15, 2019
# Make Exonerate executables for all slit fasta files (cds2genome). 
# To run: python Exonerate_multiFile.py -g path_to_genome -a path_to_anotation -f Dir_with_multiple_fasta



import os
import subprocess
import datetime
import glob
import argparse, os, pysam, sys


def parse_args():
    parser = argparse.ArgumentParser(description='Create executables for exonerate Exonerate (cds2genome)')
    parser.add_argument('-f', '--fastaDir', metavar='STR', help='Directory with all split CDS fasta files', type=str, default='Fasta_CDS/')
    parser.add_argument('-a', '--annotation', metavar='STR', help='CDS annotation file (.tab)', type=str, default='')
    parser.add_argument('-m', '--moduleLoad', metavar='STR', help='Commands to load exonerate', type=str, default="module load python/2.7-anaconda-4.4 && source activate Cori_new && ")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-g', "--genome", help="Path to genome", type=str)
    args = parser.parse_args()
    

    if (args.genome.endswith('.fa') | args.genome.endswith('.fasta')):
        print "Input genome:", args.genome
    else:
        print "Please provide a genome file"
        sys.exit(2) 
    return(args)
    
    
    
def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
        print "made new dir: ", file_path
    return()


def make_prot2genome_exonerate_executable(path_with_ALL_single_fasta_files, genome_dir, anot_file, MODULE_LOAD):
    test ="False"
    today_date = datetime.date.today()
    
    ## Parisng user input
    path_with_ALL_single_fasta_files = ''.join((path_with_ALL_single_fasta_files, "*fa"))
    fasta_single_files =  sorted(list(glob.glob(path_with_ALL_single_fasta_files)))
    genome_name        =  genome_dir.split("/")[-1].split(".fa")[0]
    exonerate_dir_out  =  ''.join(("EXONERATE_",genome_name,"/"))
    Exec_out_file      =  ''.join((genome_name,"_",str(today_date),"_exonerate.sh"))
    
    
    ex_com_1  = "exonerate --model coding2genome"
    ex_com_2 = "--score 150 --percent 20 --softmaskquery no --softmasktarget no --showcigar yes  --showvulgar yes --minintron 20 --maxintron 1000000 --geneseed 50 --showtargetgff yes --showquerygff yes --ryo"
    ex_ryo = '"[%ps]\\t%S%C" '

    print "Working directory:\t%s"    %os.getcwd()
    print "Fasta sequences:\t%s"      %path_with_ALL_single_fasta_files
    n_seqs = len(fasta_single_files)
    if len(anot_file)>0:
        print "CDS annotation file:\t%s" %anot_file
    print "Genome path:\t\t%s"     %genome_dir
    print "# of fasta files\t%s " %n_seqs
    print "Generated executable:\t%s" %Exec_out_file
    print "Output directory\t%s"      %exonerate_dir_out
    
    if not os.path.exists(exonerate_dir_out):
        os.makedirs(exonerate_dir_out)
        print "Created directory for:\t\t\n", exonerate_dir_out
    

    
    F_out = open(Exec_out_file,"w") 

    
    for single_fasta in fasta_single_files:
        sequence_name      = single_fasta.split("/")[-1].split(".fa")[0]
        seq_query          = single_fasta
        output_sinle_fasta = ''.join((">",exonerate_dir_out,sequence_name,".exonerate"))
        touch_file = ''.join((exonerate_dir_out,sequence_name,'.done'))
        touch_final        = ''.join((' && touch ', touch_file))
        if not os.path.exists(touch_file):
            if len(anot_file) > 0:
                annot_command      = ' '.join(("--annotation",anot_file))
                exonerate_commands = ' '.join((MODULE_LOAD, ex_com_1, seq_query, genome_dir, ex_com_2, ex_ryo, annot_command, output_sinle_fasta,touch_final,"\n"))
            else:
                exonerate_commands = ' '.join((MODULE_LOAD, ex_com_1, seq_query, genome_dir, ex_com_2, ex_ryo, output_sinle_fasta,touch_final,"\n"))
            F_out.write(exonerate_commands)
    F_out.close()
    
    ####HERE"
    return(Exec_out_file)


def run_executable(Exec_out_file):
    n_commands = file_len(Exec_out_file)
    qbatch_folder =  Exec_out_file.split(".sh")[0]
    
    print "Your executive file has #%s commands " %n_commands
    if n_commands > 10 : 
        #num_rounds    = int(round(n_commands/10))
        num_rounds = 5
        time          = "23:50:00"
        exeq_run      =  ' '.join(('qbatch submit -p 1 -R 3 -t',time, "-S",str(num_rounds), "-n", qbatch_folder, Exec_out_file, qbatch_folder ))
    else:
        time          = "00:30:00"
        exeq_run      = ' '.join(('qbatch submit -p 1 -R 3 -t', time, "-n", qbatch_folder, Exec_out_file, qbatch_folder ))
    return(exeq_run)

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def main():
    args = parse_args()
    MODULE_LOAD =  args.moduleLoad
#    MODULE_LOAD = "module load python/2.7-anaconda-4.4 && source activate Cori_new && "
    Exec_out_file = make_prot2genome_exonerate_executable(args.fastaDir, args.genome, args.annotation, MODULE_LOAD)
    exeq_run      = run_executable(Exec_out_file)   
    print exeq_run



if __name__=='__main__':
    main()
