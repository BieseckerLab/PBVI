# PBVI

PBVI
Position-Based Variant Identification (PBVI) is a method used to identify low-level variants by directly modeling sequencing error from a control dataset.

This repository contains model files and scripts that can be used to detect low-level variants in 11 overgrowth syndrome genes. It also contains custom scripts used to assess the sensitivity of this approach.

The PBVI scripts are written in Python 2 (version 2.7) and are not compatible with Python 3.

## Overview of PBVI

The PBVI method involves two steps. First, a control population dataset is used to generate a model that can be used for variant calling. Second, variants in input BAM files are called using Fisher’s exact test.

### 1. Model building

The model can be built by modifying variables (see below) in the generate_OGPB_model.py script. The following data are required for this step:

* Control population BAM files: the BAM files with which to build the model.
* BED file: the genomic locations to examine.


> bamfiles = glob.glob('/BAM_FILE_DIRECTORY/*.bam') #CHANGE THIS TO THE DIRECTORY WHERE BAM FILES ARE, this will look for \*bam files
> bedfile = 'OGPB_model_location.bed' #CHANGE THIS TO INPUT BEDFILE

The generate_OGPB_model.py script will create a data frame of aggregated counts and output this to a MODEL.txt file.

### 2. Variant calling

After building the model, the PBVI.py script can be modified (see below) to call variants in input BAM files. The following data are required for this step:

*Input BAM files: the BAM files in which to call variants.
*BED file: the BED file used to build the model.
*Reference genome file: a FASTA reference file.
*Model text file: the MODEL.txt file output from the generate_OGPB_model.py script.


> bamfiles = glob.glob('/BAM_FILE_DIRECTORY/*.bam') #CHANGE THIS TO THE DIRECTORY WHERE BAM FILES ARE, this will look for *bam files
> bedfile = 'OGPB_model_location.bed' #CHANGE THIS TO INPUT BEDFILE
> fastafile = 'hg19.fa' #CHANGE THIS TO FASTA REFERENCE FILE
> posmodelfile = 'OGPB_model.txt' #CHANGE THIS TO MODEL FILE

For each input BAM file, the PBVI.py script will output a text file containing the called variants.

The generate_OGPB_model.py and PBVI.py scripts are available in PBVO-OGPB directory.
