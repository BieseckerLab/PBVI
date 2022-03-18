#!/usr/bin/env python

#Attribution
#This software was developed by Jeffrey N. Dudley at the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH). Please include proper attribution of the NHGRI as the developer of this program and include a link to the following [https://github.com/BieseckerLab/PBVI] in all publications and other public disclosures that reference the program and/or include data or research results that were generated utilizing the program.

#Public Domain Notice
#This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose. No indemnification is intended or provided by the US government.

import os
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict
from functools import partial
import pysam
import time
pd.set_option('display.max_columns', None)

def get_alt_pvals(df, posmodel, CountThresh=0, alternative='two-sided', method='stouffer'):
    posmodel_df = pd.read_csv(posmodel, sep='\t', index_col=[0, 1])
    df.set_index(['Chromosome', 'Position'], inplace=True)
    df = df[df['AltCount'] > CountThresh]
    df = df.join(posmodel_df, lsuffix='_var', rsuffix='_posmodel')
    
    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_choices = [df['A_var'], df['C_var'], df['G_var'], df['T_var']]
    model_choices = [df['A_posmodel'], df['C_posmodel'], df['G_posmodel'], df['T_posmodel']]
    df['TestPosRef'] = np.select(ref_conditions, var_choices)
    df['TestPosAlt'] = np.select(alt_conditions, var_choices)
    df['ModelPosRef'] = np.select(ref_conditions, model_choices)
    df['ModelPosAlt'] = np.select(alt_conditions, model_choices)
    
    df.dropna(subset=['ModelPosRef'], inplace=True)
    df['PosPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosAlt, r.ModelPosAlt], [r.TestPosRef, r.ModelPosRef]], alternative=alternative), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosPval'] = df['PosPval'].astype(float)

    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_for_choices = [df['A_for_var'], df['C_for_var'], df['G_for_var'], df['T_for_var']]
    var_rev_choices = [df['A_rev_var'], df['C_rev_var'], df['G_rev_var'], df['T_rev_var']]
    model_for_choices = [df['A_for_posmodel'], df['C_for_posmodel'], df['G_for_posmodel'], df['T_for_posmodel']]
    model_rev_choices = [df['A_rev_posmodel'], df['C_rev_posmodel'], df['G_rev_posmodel'], df['T_rev_posmodel']]
    df['TestPosForRef'] = np.select(ref_conditions, var_for_choices)
    df['TestPosForAlt'] = np.select(alt_conditions, var_for_choices)
    df['ModelPosForRef'] = np.select(ref_conditions, model_for_choices)
    df['ModelPosForAlt'] = np.select(alt_conditions, model_for_choices)
    df['PosForPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosForAlt, r.ModelPosForAlt], [r.TestPosForRef, r.ModelPosForRef]], alternative=alternative),axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosForPval'] = df['PosForPval'].astype(float)

    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_for_choices = [df['A_for_var'], df['C_for_var'], df['G_for_var'], df['T_for_var']]
    var_rev_choices = [df['A_rev_var'], df['C_rev_var'], df['G_rev_var'], df['T_rev_var']]
    model_for_choices = [df['A_for_posmodel'], df['C_for_posmodel'], df['G_for_posmodel'], df['T_for_posmodel']]
    model_rev_choices = [df['A_rev_posmodel'], df['C_rev_posmodel'], df['G_rev_posmodel'], df['T_rev_posmodel']]
    df['TestPosRevRef'] = np.select(ref_conditions, var_rev_choices)
    df['TestPosRevAlt'] = np.select(alt_conditions, var_rev_choices)
    df['ModelPosRevRef'] = np.select(ref_conditions, model_rev_choices)
    df['ModelPosRevAlt'] = np.select(alt_conditions, model_rev_choices)
    df['PosRevPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosRevAlt, r.ModelPosRevAlt], [r.TestPosRevRef, r.ModelPosRevRef]], alternative=alternative), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosRevPval'] = df['PosRevPval'].astype(float)

    df['PosCombinePval'] = df.apply(lambda r: stats.combine_pvalues([r.PosForPval, r.PosRevPval], method=method), axis=1).str[1]

    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    ref_for_choices = [(df['A_for_var']), (df['C_for_var']), (df['G_for_var']), (df['T_for_var'])]
    df['TestPosForRef'] = np.select(ref_conditions, ref_for_choices)
    ref_rev_choices = [(df['A_rev_var']), (df['C_rev_var']), (df['G_rev_var']), (df['T_rev_var'])]
    df['TestPosRevRef'] = np.select(ref_conditions, ref_rev_choices)
    df['SBPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosForAlt, r.TestPosRevAlt], [r.TestPosForRef, r.TestPosRevRef]]), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    
    return df

def get_alt_stats(df, Depth='Depth', Ref='Ref'):
    df[Ref] = df[Ref].str.upper()
    df['freqA'] = df['A'] / df[Depth]
    df['freqC'] = df['C'] / df[Depth]
    df['freqG'] = df['G'] / df[Depth]
    df['freqT'] = df['T'] / df[Depth]
    df['freqAltA'] = np.where(df[Ref] == 'A', 0, df['freqA'])
    df['freqAltC'] = np.where(df[Ref] == 'C', 0, df['freqC'])
    df['freqAltG'] = np.where(df[Ref] == 'G', 0, df['freqG'])
    df['freqAltT'] = np.where(df[Ref] == 'T', 0, df['freqT'])
    df['AltACount'] = np.where(df[Ref] == 'A', 0, df['A'])
    df['AltCCount'] = np.where(df[Ref] == 'C', 0, df['C'])
    df['AltGCount'] = np.where(df[Ref] == 'G', 0, df['G'])
    df['AltTCount'] = np.where(df[Ref] == 'T', 0, df['T'])
    df['freqAlt'] = df[['freqAltA', 'freqAltC', 'freqAltG', 'freqAltT']].max(axis=1)
    df['Alt'] = df[['freqAltA', 'freqAltC', 'freqAltG', 'freqAltT']].idxmax(axis=1).str[7:]
    conditions = [(df['Alt'] == 'A'), (df['Alt'] == 'C'), (df['Alt'] == 'G'), (df['Alt'] == 'T')]
    choices = [df['A'], df['C'], df['G'], df['T']]
    df['AltCount'] = np.select(conditions, choices)
    df['AltCount'] = np.where(df['freqAlt'] == 0, 0, df['AltCount'])
    return df

def get_basecounts(bam, bed, fasta):
    def check_read(read, nodup=True, both=True, forward=False):
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

    def check_read_mapqual(read, nodup=True, both=True, forward=False):
        if not read.is_unmapped and not read.is_qcfail and not read.is_secondary and read.mapping_quality >= 20:
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

    def split_chr_df(counts):
        chr_list = [{k: v} for (k, v) in counts.iteritems()]
        dfs = []
        for chrom in chr_list:
            name = list(chrom.iterkeys())[0]
            vals = chrom[name]
            df = pd.DataFrame.from_dict(data=vals, orient='index', columns=['A', 'C', 'G', 'T'])
            df['Chromosome'] = name
            df['Depth'] = df['A'] + df['C'] + df['G'] + df['T']
            df.index = pd.MultiIndex.from_tuples(df.index, names=['Position', 'Ref'])
            df.reset_index(drop=False, inplace=True)
            dfs.append(df)
        df = pd.concat(dfs)
        df.set_index(['Chromosome', 'Position', 'Ref'], inplace=True)
        return df

    bed_list = []
    for_count_dict = defaultdict(lambda: defaultdict(partial(np.zeros, 4, dtype=int)))
    rev_count_dict = defaultdict(lambda: defaultdict(partial(np.zeros, 4, dtype=int)))
    b = 0
    fa = pysam.FastaFile(fasta)
    with open(bed) as bedfile:
        print('now counting %s' % os.path.basename(bam))
        samfile = pysam.AlignmentFile(bam, 'rb')
        for line in bedfile:
            chromo = line.split('\t', 10)[0]
            start = int(line.split('\t', 10)[1])
            end = int(line.split('\t', 10)[2])
            if not (chromo, start, end) in bed_list:
                bed_list.append((chromo, start, end))
                x = 0
                y = 0
                ref_bases = fa.fetch(chromo, start, end)
                region = np.arange(start+1, end + 1, 1)
                for nuc in samfile.count_coverage(chromo, start, end, quality_threshold=0, read_callback=(lambda z: check_read_mapqual(z, both=False, forward=True))):
                    for idx, count in enumerate(nuc):
                        for_count_dict[chromo][(region[idx], ref_bases[idx].upper())][x] = count
                    x += 1
                for nuc in samfile.count_coverage(chromo, start, end, quality_threshold=0, read_callback=(lambda z: check_read_mapqual(z, both=False, forward=False))):
                    for idx, count in enumerate(nuc):
                        rev_count_dict[chromo][(region[idx], ref_bases[idx].upper())][y] = count
                    y += 1
    df = split_chr_df(for_count_dict).join(split_chr_df(rev_count_dict), lsuffix='_for', rsuffix='_rev')
    df['A'] = df['A_for'] + df['A_rev']
    df['C'] = df['C_for'] + df['C_rev']
    df['G'] = df['G_for'] + df['G_rev']
    df['T'] = df['T_for'] + df['T_rev']
    df['Depth'] = df['Depth_for'] + df['Depth_rev']
    df.reset_index(inplace=True)
    return df


bamfiles = glob.glob('/BAM_FILE_DIRECTORY/*.bam') #CHANGE THIS TO THE DIRECTORY WHERE BAM FILES ARE, this will look for *bam files

bedfile = 'OGPB_model_location.bed' #CHANGE THIS TO INPUT BEDFILE

fastafile = 'hg19.fa' #CHANGE THIS TO FASTA REFERENCE FILE

posmodelfile = 'OGPB_model.txt' #CHANGE THIS TO MODEL FILE


var_dfs = []
alternative = 'greater'
method = 'stouffer'
for bamfile in bamfiles:
    print(bamfile)
    df = get_basecounts(bamfile, bedfile, fastafile)

    df.sort_values(by=['Position'], inplace=True, ascending=True)

    df.to_csv(os.path.basename(bamfile) + '_PBVI_output_unannotated.txt', sep='\t', index=True)
    
    alt_df = get_alt_stats(df)
    var_df = get_alt_pvals(df=alt_df, posmodel=posmodelfile, CountThresh=0, alternative=alternative, method=method)
    ID = os.path.basename(bamfile).split('.', 10)[3] # Might need to modify this depending on BAM file names
    var_df['ID'] = ID
    
    var_df.to_csv(os.path.basename(bamfile)+'_PBVI_output_annotated.txt', sep='\t', index=True)


