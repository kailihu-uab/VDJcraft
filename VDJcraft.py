import argparse
import re
#from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import sys
import shutil
import pandas as pd
import script
import dscript
import high_confident
import D_extract
import error_correction
import io
import cdr3new
import cdr1_cdr2
import tempfile
import subprocess
import os
import math
#need to install 'pip install biopython' to run import blast
#f = open('output.txt', 'w')
#parser = argparse.ArgumentParser(description="print Hello")

parser = argparse.ArgumentParser(description="blastn")

#set up the blastn command line

parser.add_argument('-I','--fastq', type=str, help="Enter your input file")
parser.add_argument('-ontI','--ontfq', type=str, help="Enter your input file")
parser.add_argument('-R','--database', type=str, help="Enter your database",required=True)
parser.add_argument('-O', '--outs', type=str, help="Set output folder",required=True)

parser.add_argument('-of','--outfile',type=str, help="output result file",required=True)

parser.add_argument('-p', '--VJIdent', type=int, default = 90, help="Enter preferred identifation percentage for VJ genes")
parser.add_argument('-pd', '--DIdent', type=int, default = 95, help="Enter preferred identifation percentage for D genes")
parser.add_argument('-s', '--score', type=int, default = 300, help="Enter preferred matching score for V genes")


args = parser.parse_args()

output_folder = args.outs

def generate_output_file(file_name):
    # Your logic to generate the output file goes here
  #  content = f"This is the content of {file_name}"

    # Specify the folder where you want to save the files
 #   output_folder = args.outs

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Construct the full path to the output file
   # file_path = os.path.join(output_folder, file_name)
    return os.path.join(output_folder, file_name)

#def run_samtools_fasta(input_file, output_file, out_table):
    # Run minimap2 to align reads to a reference genome
#    output_file="out.fa"

in_file = args.fastq
#ont_file=args.ontfq
temp = "temp_folder"
temp_dir = generate_output_file(temp)
os.makedirs(temp_dir, exist_ok=True)
 #   try:
minimap2_output_file = os.path.join(temp_dir, 'minimap2_output.sam')
minimap2_cmd = ['minimap2', '-t 8', '-ax', 'splice:hq', '-uf', '/data/scratch/klhu0502/raw_data/GRCh38_full_analysis_set_plus_decoy_hla.fa', in_file]
#minimap2_cmd = ['minimap2', '-t 8', '-ax', 'splice','/data/scratch/klhu0502/raw_data/GRCh38_full_analysis_set_plus_decoy_hla.fa', ont_file]
with open(minimap2_output_file, 'w') as minimap2_output:
    minimap2_proc = subprocess.run(minimap2_cmd, stdout=minimap2_output)
      #      minimap2_proc = subprocess.run(minimap2_cmd, stdout=subprocess.PIPE)

    # Run samtools to convert the minimap2 output to BAM format
samtools_bam_file = os.path.join(temp_dir, 'samtools_output.bam')
   #     samtools_cmd = ['samtools', 'view', '-b', '-']
samtools_cmd = ['samtools', 'view', '-b', minimap2_output_file]
with open(samtools_bam_file, 'wb') as samtools_bam_output:
    samtools_proc = subprocess.run(samtools_cmd, stdout=samtools_bam_output)
      #  samtools_proc = subprocess.run(samtools_cmd, input=minimap2_proc.stdout, stdout=subprocess.PIPE)

sorted_bam_file = os.path.join(temp_dir, 'sorted_output.bam')
samtools_cmd_sort = ['samtools', 'sort', samtools_bam_file]
with open(sorted_bam_file, 'wb') as sorted_bam_output:
    samtools_proc_sort = subprocess.run(samtools_cmd_sort, stdout=sorted_bam_output)

samtools_index_file = sorted_bam_file + '.bai'
samtools_cmd_index = ['samtools', 'index', sorted_bam_file]
samtools_proc_index = subprocess.run(samtools_cmd_index)

index_file_destination = os.path.join(temp_dir, os.path.basename(samtools_index_file))
shutil.move(samtools_index_file, index_file_destination)

samtools_vdj_bam = os.path.join(temp_dir, 'samtools_vdj.bam')
samtools_cmd3 = ['samtools', 'view', '-b', '-L', '/data/user/klhu0502/work/longscript/vdjtool/VDJC_genome_CDSuniq.bed', sorted_bam_file]
with open(samtools_vdj_bam, 'wb') as vdj_bam:
    samtools_proc_vdj = subprocess.run(samtools_cmd3, stdout=vdj_bam)

samtools_cmd_view = ['samtools', 'view', samtools_vdj_bam]
samtools_proc_view = subprocess.run(samtools_cmd_view, stdout=subprocess.PIPE)
cut_cmd = ['cut', '-f', '1,3,4']
cut_proc = subprocess.run(cut_cmd, input=samtools_proc_view.stdout, stdout=subprocess.PIPE)
columns = ['Query', 'chr', 'position']
pos_df = pd.read_csv(io.StringIO(cut_proc.stdout.decode()), sep='\t', header=None, names=columns)

#fastq = args.fastq
#in_file = args.fastq
table ='table.csv'
tablefile=generate_output_file(table)
with open (tablefile, 'w') as c:
    pos_df.to_csv(c, sep='\t', index=False, na_rep='NA')

output_file ='vdj_candidate.fa'
file_path=generate_output_file(output_file)
samtools_cmd5 = ['samtools', 'fasta', '-F 4', samtools_vdj_bam]
with open (file_path, 'w') as k:
    subprocess.run(samtools_cmd5, stdout=k)

k.close()
c.close()
vdj_bam.close()
sorted_bam_output.close()
minimap2_output.close()

shutil.rmtree(temp_dir)

blastout = "blastout.txt"
file_blast=generate_output_file(blastout)
if os.path.exists(file_blast):
    os.remove(file_blast)

#dseq_file ='/data/user/klhu0502/work/longscript/dseq.fa'

blastn_cmd = ['blastn','-query',file_path,'-db',args.database,'-out',file_blast]
subprocess.run(blastn_cmd)

#blastn_cline = NcbiblastnCommandline(query=args.fasta, db=args.database, out=file_blast)
#blastn_d_cline = NcbiblastnCommandline(query=output_file, db=args.database, word_size = 7, evalue = 1000, out=blastoutd)
#stdout, stderr = blastn_cline()
#stdout1, stderr1 = blastn_d_cline()


#if os.path.exists('-.bai'):
#    os.remove('-.bai')

output_file = "raw_gene.txt"
raw=generate_output_file(output_file)
if os.path.exists(raw):
    os.remove(raw)

with open(file_blast, 'r') as f:
    script.blastfmt(f,raw)

f.close()

addline='addline.csv'
add_line=generate_output_file(addline)
vdjout="vdjout.csv"
vdjc_out=generate_output_file(vdjout)
vjfile="vjtemp.csv"
vjtemp=generate_output_file(vjfile)
vgene="vgene.csv"
vgene_temp=generate_output_file(vgene)
#with open(raw, 'r') as g:
 #   lines = g.readlines()
high_confident.top_extract(raw,add_line,vdjc_out,vjtemp,vgene_temp)
#g.close()

dseq_file ='dgene_candidate.fa'
file_dseq =generate_output_file(dseq_file)

blastoutd = "blastoutd.txt"
file_blastd=generate_output_file(blastoutd)
if os.path.exists(file_blastd):
   os.remove(file_blastd)

#/data/user/klhu0502/work/longscript/dgenedb/IMGTD.fa
with open(vjtemp,'r') as vj:
    bed_vj=pd.read_csv(vjtemp, sep='\t', header=0)
    bed_vj['dstart'] = bed_vj.apply(lambda row: row['v_queryend']-20 if row['strand'] == 'Plus/Plus' else row['j_queryend']-20, axis=1)
    bed_vj['dend'] = bed_vj.apply(lambda row: row['j_querystart']+20 if row['strand'] == 'Plus/Plus' else row['v_querystart']+20, axis=1)
    extracted_sequences = D_extract.extract_sequences(file_path, bed_vj, file_dseq)
    blastn_d = ['blastn','-query',file_dseq,'-db','/data/user/klhu0502/work/longscript/dgenedb/IMGTD.fa','-word_size','7','-evalue','0.05','-out',file_blastd]
    subprocess.run(blastn_d)

vj.close()

output_filename = "raw_dgene.txt"
rawd=generate_output_file(output_filename)

with open(file_blastd, 'r') as d:
    dscript.blastdgene(d,rawd)

d.close()

dout="top_dgene.csv"
d_out=generate_output_file(dout)

#with open(rawd, 'r') as dfile:

  #  lines = dfile.readlines()

D_extract.top_extract(rawd,d_out)

#dfile.close()

orig_vdj = "vdjc_original.csv"
mergevdjc= generate_output_file(orig_vdj)

with open(vdjc_out,'r') as vjc, open(d_out,'r') as d:

    df_vjc = pd.read_csv(vjc, sep='\t', header=0)

    df_d = pd.read_csv(d, sep='\t', header=0)

    merge_vdjc = pd.merge(df_vjc, df_d, on = 'Query', how='outer')

    merge_vdjc = merge_vdjc[['v_gene','d_gene','j_gene','c_gene','v_seq','d_seq','j_seq','Query']]

#    sort_vdjc = merge_vdjc.iloc[1:].sort_values(by=merge_vdjc.columns[0])

    merge_vdjc.to_csv(mergevdjc, sep='\t', index=False, na_rep='OFM')

vjc.close()
d.close()

corr_vdjc = "vdjc_corrected.csv"
correct_vdjc = generate_output_file(corr_vdjc)
with open(mergevdjc,'r') as vdjc, open(tablefile,'r') as T:
    df_vdjc = pd.read_csv(vdjc, sep='\t', header=0)

    df_pos = pd.read_csv(T, sep='\t', header=0)

    pos_vdjc = pd.merge(df_pos, df_vdjc, on = 'Query', how='inner')

    corrected_df = error_correction.error_correction(pos_vdjc)

    sorted_df = corrected_df.iloc[1:].sort_values(by=corrected_df.columns[0])

    sorted_df.to_csv(correct_vdjc, sep='\t', index=False)

T.close()


cdr3 = "cdr3.csv"
cdr3_temp = generate_output_file(cdr3)
with open(vjtemp,'r') as a:
    df_vjtemp = pd.read_csv(a, sep='\t', header=0)
    cdr3 = cdr3new.cdr3_extract(df_vjtemp,file_path)
    cdr3.to_csv(cdr3_temp, sep='\t', index=False)

a.close()

cdr1 = "cdr1.csv"
cdr1_temp = generate_output_file(cdr1)
cdr2 = "cdr2.csv"
cdr2_temp = generate_output_file(cdr2)
with open(vgene_temp,'r') as cd:
    df_v = pd.read_csv(cd, sep='\t', header=0)
    cdr1 = cdr1_cdr2.extract_cdr1_sequence(df_v)
    cdr1.to_csv(cdr1_temp, sep='\t', index=False)
    df2=pd.read_csv('/data/user/klhu0502/work/longscript/vdjtool/CDR2.consensus.csv', sep='\t',header=0)
    #*mergecdr2 = pd.merge(merge_vdjc, df2, on = 'v_gene')
    df_v['v_gene'] = df_v['v_gene'].str.replace(r'_.*', '', regex=True)
    df_cdr2 = pd.merge(df_v, df2, on = 'v_gene', how='inner')
    cdr2 = cdr1_cdr2.cdr2_extract(df_cdr2, 'v_gene','strand','v_seq','cdr2seq','cdr2_rev','Query')
    cdr2.to_csv(cdr2_temp, sep='\t', index=False)

cd.close()



