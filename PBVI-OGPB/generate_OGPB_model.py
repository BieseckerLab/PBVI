#!/usr/bin/env python

#Attribution
#This software was developed by Jeffrey N. Dudley at the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH). Please include proper attribution of the NHGRI as the developer of this program and include a link to the following [https://github.com/BieseckerLab/PBVI] in all publications and other public disclosures that reference the program and/or include data or research results that were generated utilizing the program.

#Public Domain Notice
#This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose. No indemnification is intended or provided by the US government.


import glob
import numpy as np
import pysam
from collections import defaultdict
from functools import partial
import os
import pandas as pd

bamfiles = glob.glob('/data/NHGRIgeno/celine/CLIA_*.hg38.markduplicate.sort.bam.recal.bam') #CHANGE THIS TO THE DIRECTORY WHERE BAM FILES ARE, this will look for *bam files
bedfile = 'OGPB_model_location.bed' #CHANGE THIS TO INPUT BEDFILE


def get_chr_df(chr_list, read='one'):
    dfs = []
    for chrom in chr_list:
        name = list(chrom.iterkeys())[0]
        vals = chrom[name]
        df = pd.DataFrame.from_dict(data=vals, orient='index', columns=['A', 'C', 'G', 'T', 'Depth'])
        df['chr'] = name
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def split_chr_df(counts):
    chr_list = [{k: v} for (k, v) in counts.iteritems()]
    dfs = []
    for chrom in chr_list:
        name = list(chrom.iterkeys())[0]
        vals = chrom[name]
        df = pd.DataFrame.from_dict(data=vals, orient='index', columns=['A', 'C', 'G', 'T', 'Depth'])
        df['Chromosome'] = name
        df.index.names = ['Position']
        df.reset_index(drop=False, inplace=True)
        dfs.append(df)
    df = pd.concat(dfs)
    df.set_index(['Chromosome', 'Position'], inplace=True)
    return df


def check_read(read, nodup=True, both=True, forward=True):
    if not read.is_unmapped and not read.is_qcfail and not read.is_secondary:
        if read.is_duplicate and nodup:
            return False
        elif both:
            return True
        elif forward and not read.is_reverse:
            return True
        elif not forward and read.is_reverse:
            return True
    else:
        return False


def model_counts(bams, bed):
    for_count_dict = defaultdict(lambda: defaultdict(partial(np.zeros, 5, dtype=int)))
    rev_count_dict = defaultdict(lambda: defaultdict(partial(np.zeros, 5, dtype=int)))
    b = 0
    for bam in bams:
        b += 1
        if b < 1000:
            with open(bed) as bedfile:
                bed_list = []
                print('now working on %s which is the %d bam' % (os.path.basename(bam), b))
                samfile = pysam.AlignmentFile(bam, 'rb')
                for line in bedfile:
                    chromo = line.split('\t', 10)[0]
                    start = int(line.split('\t', 10)[1])
                    end = int(line.split('\t', 10)[2])
                    if not (chromo, start, end) in bed_list:
                        bed_list.append((chromo, start, end))
                        x = 0
                        y = 0
                        region = np.arange(start, end + 1, 1)
                        depth_dict = defaultdict(int)
                        for column in samfile.pileup(chromo, start - 1, end):
                            rpos = column.reference_pos + 1
                            if rpos in region:
                                depth_dict[rpos] = column.nsegments
                        for nuc in samfile.count_coverage(chromo, start - 1, end, quality_threshold=0, read_callback=(
                        lambda z: check_read(z, both=False, forward=True))):
                            for idx, count in enumerate(nuc):
                                # always add to depth counts
                                for_count_dict[chromo][region[idx]][4] += count
                                for_count_dict[chromo][region[idx]][x] += count
                            x += 1
                        for nuc in samfile.count_coverage(chromo, start - 1, end, quality_threshold=0, read_callback=(
                        lambda z: check_read(z, both=False, forward=False))):
                            for idx, count in enumerate(nuc):
                                # always add to depth counts
                                rev_count_dict[chromo][region[idx]][4] += count
                                rev_count_dict[chromo][region[idx]][y] += count
                            y += 1
    df = split_chr_df(for_count_dict).join(split_chr_df(rev_count_dict), lsuffix='_for', rsuffix='_rev')
    df['A'] = df['A_for'] + df['A_rev']
    df['C'] = df['C_for'] + df['C_rev']
    df['G'] = df['G_for'] + df['G_rev']
    df['T'] = df['T_for'] + df['T_rev']
    df['Depth'] = df['Depth_for'] + df['Depth_rev']
    return df


df = model_counts(bamfiles, bedfile)

print(df)


df.to_csv('MODEL.txt', sep='\t')
