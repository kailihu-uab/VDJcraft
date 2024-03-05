# VDJcraft
 comprehensively characterze immune repertoire

# V(D)Jcraft

A V(D)J identifer for long-read transcriptome data.

Author: Kaili Hu

Email: klhu0502@uab.edu

Draft date: Feb.12, 2024

## Quick Start
```sh
git clone .git
cd VDJcraft/
./VDJcraft -h

To simplify the environment setup process, Anaconda2 (https://www.anaconda.com/) is recommended. To create an environment with conda:
conda create --name vdjcraft python=3
conda activate vdjcraft
conda install bioconda::minimap2
conda install bioconda::samtools
conda install bioconda::blast

pip3 install biopython
pip install python-Levenshtein

vdjcraft.py -h

# quick using VDJcraft with fastq
vdjcraft -I input.fastq -R IMGT.fa -O output_folder -of output

# VDJcraft discovery with custom parameter
vdjcraft -I input.fastq -R IMGT.fa -O output_folder -of output
```

## Description

V(D)Jcraft is a tool for VDJC identification with long-read transcriptome sequencing data. The input should be sequencing reads file: fastq or fa (PacBio Iso-Seq, Nanopore, or mixed platform). The output is a list of confident VDJC gene list and their
 CDR nt sequences and amino acid sequences. By default, V(D)Jcraft uses Human GRCh38 and IMGT VDJC database.
#When using custom reference, make sure the chromosome name in BAM, GTF, and reference_genome.fa are identical. By default, FusionSeeker only considers gene with valid "gene_name" in GTF and skips the remaining genes, unless --geneid is set.<br />
This program was tested on a x86_64 Linux system with a 128GB physical memory.


## Depencency

Dependencies for FusionSeeker:

* python3
* minimap2  (tested with version 2.24)
* samtools  (tested with version 1.9)
* blastn


## Installation

```
git clone https://github.com/Kaili-Hu/V(D)Jcraft.git
```
Then, please also add this directory to your PATH:
```
export PATH=$PWD/VDJcraft/:$PATH
```


To simplify the environment setup process, Anaconda2 (https://www.anaconda.com/) is recommended:
```
conda create --name vdjcraft -y
conda activate vdjcraft
conda install -c bioconda minimap2=2.24 pysam=0.17 samtools=1.9 -y

```

A test dataset is available to verify successful installation:
```
VDJcraft -I vdjcraft/testdata/test.fastq  -R vdjcraft/testdata/IMGT.fa -o test_out/ -of test_output
```
Output should be identical to confident_vdjc.txt in the testdata folder, with 30 vdjc.
(The VDJC identification on test dataset should finish within several minutes with 1 CPUs and 1GB memory.)

## General usage


```
VDJcraft [-h] -I <input.fastq>

Gene fusion caller for long-read sequencing data

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --I Fastq             Input fastq.
  --datatype DATATYPE   Input read type (isoseq, nanopore) [nanopore]
  --gtf GTF             Genome annotation file
  --R REF               Reference genome. Required for breakpoint polishing
  --human38             Use reference genome and GTF for Human GCRh38 (default)
  -o OUTPATH, --outpath OUTPATH
                        Output directory [./fusionseeker_out/]
  -of OUTPUT            Output filename

```

## Use cases
VDJcraft requires a input of sequencing reads (Fasta or Fastq format):
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
VDJcraft -I vdjcraft/testdata/test.fastq  -R vdjcraft/testdata/IMGT.fa -o test_out/ -of test_output
```
Or with custom reference genome and filteration cutoff:
```
VDJcraft -I vdjcraft/testdata/test.fastq  -R vdjcraft/testdata/IMGT.fa -o test_out/ -of test_output
```

By default, VDJcraft uses filteration with identification percentage 90%, and select the hightest score from blast output:
```
VDJcraft -I vdjcraft/testdata/test.fastq  -R vdjcraft/testdata/IMGT.fa -o test_out/ -of test_output
```



