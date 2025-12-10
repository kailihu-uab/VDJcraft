# VDJcraft

A V(D)J identifer for long-read transcriptome data.

Author: Kaili Hu

Email: klhu0502@uab.edu

Draft date: Dec.2rd, 2025

## Quick Start
```sh
git clone git clone https://github.com/kailihu-uab/VDJcraft.git
cd VDJcraft/
chmod +x vdjcraft
./VDJcraft -h

# quick using VDJcraft with fastq
./VDJcraft -I input.fastq -ref genome.fa -O output_folder -of output

# VDJcraft discovery with custom parameter
./VDJcraft -I input.fastq -ref T2T.fa -t2t -p 95 -pd 70 -O output_folder -of output

# VDJcraft discovery of potential novel events
./VDJcraft -I test.fastq -shm -ref GRCh38.p14.genome.fa -O output_folder -of output

# VDJcraft discovery with other references
./VDJcraft -I test.fastq -hs37 -ref GRCh37.genome.fa -O output_folder -of output

## Description

V(D)Jcraft is a tool for VDJC identification with long-read transcriptome sequencing data. The input should be sequencing reads file: fastq or fa (PacBio Iso-Seq, Nanopore, or mixed platform). The output is a list of confident VDJC gene list and their CDR nt sequences and amino acid sequences. By default, V(D)Jcraft uses Human GRCh38 and IMGT VDJC database.
#When using custom reference, make sure the chromosome name in BAM, and reference_genome.fa are identical. By default, VDJcraft use GRCh38 as reference, please download latest human reference from GENCODE.<br />
This program was tested on a x86_64 Linux system with a 128GB physical memory.


## Depencency

Dependencies for VDJcraft:

* python3 (tested with version 3.8.16)
* minimap2  (tested with version 2.26)
* samtools  (tested with version 1.17)
* blastn (tested with version 2.14.0+)


## Installation

```
git clone https://github.com/Kaili-Hu/V(D)Jcraft.git
```

Install conda:
To simplify the environment setup process, Anaconda3 (https://www.anaconda.com/) or Miniconda(https://repo.anaconda.com/miniconda/) is recommended.

On HPC/cluster systems if conda module is available:
```
module avail
module load anaconda3
```

After conda installed or loaded, run commands as below to load dependencies:
```
conda create --name vdjcraft python=3 -y
conda activate vdjcraft
conda install -c bioconda minimap2 samtools blast -y
pip install biopython python-Levenshtein

```

After installing dependencies, run commands as below:
```
cd VDJcraft/
chmod +x vdjcraft
./VDJcraft -h
```

Then,please also add this directory to your PATH (optional):
```
export PATH=$PWD/VDJcraft/:$PATH
```

A test dataset is available to verify successful installation:
```
./VDJcraft -I test.fastq -ref GRCh38.p14.genome.fa -O test -of vdjcraft_test
```
Output should be identical to vdjcraft_report.csv in the testdata folder, with 16 identification.
(The VDJC identification on test dataset should finish within several minutes with 1 CPUs and 1GB memory.)

## General usage


```
./VDJcraft [-h] -I <input.fastq>

V(D)J recombination caller for long-read sequencing data

optional arguments:
    -h, --help            show this help message and exit
    -v, --version         show program's version number and exit
    -I, Fastq             Input fastq.
    -bam, input bam       Input bam alignment file
    -ref, reference       Enter your human reference genome.
    -O, outfolder         Set output folder.
    -of, outfile          Output result file.
    -ont, ONT fq          Enter your input ONT file.
    -t2t, T2Tbed          If using T2T as human refence genome, has to use this which corresponds to T2T BED file.
    -hs37, GRCh37         If using GRCh37 as human refence genome, has to use this which corresponds to hg37 BED file.
    -Mm, mouse-GRCm39     If using GRCm39 as mouse refence genome, has to use this which corresponds to GRCm39 BED file.
    -p, VJIdent           Enter preferred identifation percentage for VJC genes, default = 90.
    -pd, DIdent           Enter preferred identifation percentage for D genes,default = 80.
    -s, score             Enter preferred matching score,default = 300.
    -shm, SHM             Enter to get novel events due to SHM.
    -m, SHM cutoff        Enter cutoff of matching rate for SHM,default = 85.

```

## Use cases
VDJcraft requires a input of sequencing reads (Fasta or Fastq format) or bam file:

```
# PacBio Iso-Seq
minimap2 -ax splice:hq reference.fa  isoseq.fastq | samtools sort -o isoseq.bam
samtools index isoseq.bam
# Nanopore
minimap2 -ax splice reference.fa  nanopore.fastq | samtools sort -o nanopore.bam
samtools index nanopore.bam
```

VDJcraft can be applied with built-in Human reference genome (hg38) and annotation (Ensembl v104):
```
./VDJcraft -I testdata/test.fastq  -ref human_refence.fa(eg.GRCh38.p14.genome.fa) -o output_folder -of output
```
Or with custom reference genome and filteration cutoff for VJC and D genes:
```
./VDJcraft -I vdjcraft/testdata/test.fastq -t2t -p 95 -pd 85 -ref T2T.fa -o output_folder -of output
```

VDJcraft can be applied with ONT fastq or fasta file with parameter '-ont':
```
./VDJcraft -ont testdata/ONT.fastq  -ref human_refence.fa(eg.GRCh38.p14.genome.fa) -o output_folder -of output
```

### Options of VDJcraft
#### 1. -shm, generating events caused by somatic hypermutation(SHM).
--shm is the most essential and biologically meaningful argument for novel candidate of VDJcraft. It is used to generate sequences potentially containning SHM and partially align to IMGT database.
By default, VDJcraft estimates partially aligning to IMGT cutoff by 85% default, customized by -m. If you find number of novel events is too high under default settings, you can speficy a lower -m cutoff to allow in less candida
tes:
```
./VDJcraft -I test.fastq -shm -m 80 -ref GRCh38.p14.genome.fa -O output_folder -of output
```

#### 2. -ref, reference genome
Input reference genome allows VDJcraft to align transcript sequences and refine VDJC positions of candidate reads extraction. Make sure to provide the correspondingparameter with same reference file used for read alignment.Defau
lt is GRCh38, custom reference can be: -t2t, T2T ref; -hg37, GRCh37; -Mm, GRCm39.
This option adjusts different referenct corresponding VDJC position bed files. By default, GRCh38 is applied, custom reference can be: -t2t, T2T ref; -hg37, GRCh37; -Mm, GRCm39

```
./VDJcraft -I test.fastq -hs37 -ref GRCh37.genome.fa -O output_folder -of output
```

#### 3. --R, IMGT database
Input IMGT database can be generated by users for second path of local realignment. Users are welcome to generate IMGT database fasta file to be blastn database. Make sure to index the custom self_generate_IMGT.fa by 'makeblastd
b -in self_generate_IMGT.fa -dbtype nucl -out self_generate_IMGT'.
By default, VDJcraft using latest version updated in 2025 on IMGT is provided.
```
./VDJcraft -I test.fastq -R self_generate_IMGT.fa -ref GRCh38.p14.genome.fa -O output_folder -of output
```


## Output files
The output directory includes:
```
vdjcraft_report.csv                     Final report of VDJcraft. A list of confident VDJC identifications from long-read transcriptome data, including count,v_gene,d_gene,j_gene,c_gene,cdr_sequence,aa_sequence
vdjc_original.txt                       A list of confident VDJ original calls before error correction from input fastq file. Includes gene names, gene sequences, name of VDJC-supporting reads.
vdjc_corrected.fa                       Reported confident VDJC gene calls after error correction.
cdr1.csv                                A full list of detected CDR1 identification and sequences, includes read name, related V genes, CDR1 sequences.
cdr2.csv                                A full list of detected CDR2 identification and sequences, includes read name, related V genes, CDR2 sequences.
cdr3.csv                                A full list of detected CDR3 identification and sequences, includes read name, related V genes, CDR3 sequences.
blastout.txt                            Results of candidate reads aligned with IMGT database.
raw_gene.txt                            A list of VJC gene raw signals after blast.
raw_dgene.txt                           A list of D gene raw signals after blast
novel.csv                               A list of potential novel events caused by SHM.
Other files                             Intermediate files during detection.
(temp_folder/                           Intermediate bam files during raw signal detection. Removed by default.)


