#!/usr/bin/env python

#Attribution
#This software was developed by Jeffrey N. Dudley at the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH). Please include proper attribution of the NHGRI as the developer of this program and include a link to the following [https://github.com/BieseckerLab/PBVI] in all publications and other public disclosures that reference the program and/or include data or research results that were generated utilizing the program.

#Public Domain Notice
#This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose. No indemnification is intended or provided by the US government.

import pandas as pd
import numpy as np
import pysam
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

from decimal import Decimal

pd.set_option('display.max_columns', None)
plt.style.use('ggplot')



def get_refs(df, fasta):
    fa = pysam.FastaFile(fasta)
    for index, row in df.iterrows():
        chromo = index[0]
        pos = index[1]
        ref = fa.fetch(chromo, pos - 1, pos)
        df.loc[index, 'Ref'] = ref
    return df


def get_model_ref_mean(df):
    conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    choices = [df['A'], df['C'], df['G'], df['T']]
    df['RefCount'] = np.select(conditions, choices)
    print('the mean ref allele count is %s' % str(np.mean(df['RefCount'])))
    refmean = int(np.mean(df['RefCount']))
    return refmean


def get_min_counts(amin, amax, astep, dmin, dmax, dstep, modelrefcount, bon):
    count_dict = defaultdict(lambda: defaultdict(float))
    strand_count_dict = defaultdict(lambda: defaultdict(float))
    allele_range = np.arange(amin, amax, astep)
    # allele_range
    depth_range = np.arange(dmin, dmax, dstep)
    depth_range=[150,300,600,1200]
    print ("in get min_counts function: allele range, depth range")
    print (allele_range)
    print (depth_range)
    
    for depth in depth_range:
        for modelcount in allele_range:
            y = 0
            thresh = False
            while not thresh:
                print ("in loop: depth and model count", depth, modelcount)
                ref_depth = int(round(depth / 2.))
                for_count = int(round(y / 2.))
                rev_count = (y - for_count)
                for_table = np.array([[for_count, int(round(modelcount / 2.))], [(ref_depth - for_count), int((round(modelrefcount / 2.)))]])
                for_pval = stats.fisher_exact(for_table, alternative='greater')[1]
                rev_table = np.array([[rev_count, int(round(modelcount / 2.))], [(ref_depth - rev_count),int((round(modelrefcount / 2.)))]])
                rev_pval = stats.fisher_exact(rev_table, alternative='greater')[1]
                strand_pval = stats.combine_pvalues([for_pval, rev_pval], method='stouffer')[1]
                print ("forward table", for_table)
                print ("reverse table", rev_table)
                if strand_pval < (0.05 / bon):
                    thresh = True
                    vaf = y/float(depth)
                    strand_count_dict[depth][modelcount] = vaf
                else:
                    y += 1

    print("strand_count", strand_count_dict)
    return count_dict, strand_count_dict



def plot_count_thresh_dict(vaf_dict, strand_dict, name):
    print ("in plot_count_thresh_dict")
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    print ("plot made")
    print(color_cycle, type(color_cycle))
    plt.figure(figsize=(8, 7))
    ax1 = plt.subplot(111)
  
    depths = sorted([d for d in strand_dict])
    print ("depths sorted")
    cs = 0
    for depth in depths:
        print ("depth", depth)
        xs = [x for x in strand_dict[depth]]
        print (xs)
        print (strand_dict[depth])
        print ("printing items", strand_dict[depth].items())
        #ys = [v for (k,v) in strand_dict[depth].iteritems()]
        ys = [v for (k,v) in strand_dict[depth].items()] ##corrected CH
        print (ys)
        #with line style
        # plt.plot(xs, ys, linewidth=3, color=color_cycle[cs], linestyle=':')
        plt.plot(xs, ys, linewidth=3, color=color_cycle[cs], label=depth)
        cs += 1
    print ("plot create")
    ax1.set_ylabel('Test minimum VAF', color='black', fontname='Arial', fontsize=16)
    ax1.set_xlabel('Model alternate nucleotide counts', color='black', fontname='Arial', fontsize=16)
    ax1.tick_params(axis='x', colors='black')
    ax1.tick_params(axis='y', colors='black')

    plt.legend(prop={'family': 'Arial'}, loc=4)
    plt.savefig(name, dpi=600)
    plt.show()



modelfile='/path/OGPB_model.txt' ##CHANGE THIS: full path to model file
fastafile='/path/hg19.fasta' ##CHANGE THIS: full path to a reference fasta file
df = pd.read_csv(modelfile, sep='\t')

df.set_index(['Chromosome', 'Position'], inplace=True, drop=True)
df = get_refs(df=df, fasta=fastafile)

print(df)

def get_alt_dist(df):
    plt.figure()
    nucs = ['A', 'C', 'G', 'T']
    alt_dict = {n: [] for n in nucs}
    print(alt_dict)
    for ref in nucs:
        ref_df = df[df['Ref'] == ref]
        alt_nucs = [x for x in nucs if x != ref]
        for alt in alt_nucs:
            print(ref, alt)
            alt_list = ref_df[alt].tolist()
            print(alt_list)
            alt_dict[alt].extend(alt_list)
    for alt in nucs:
        counts = alt_dict[alt]
        sns.distplot(counts, label=alt, hist=False, kde_kws={'linewidth': 3})
    axes = plt.gca()
    axes.set_xlim([0, 100])
    plt.legend()
    plt.savefig('alt_dists_pik3ca_nohist.png', dpi=600)
    plt.show()




refmean = get_model_ref_mean(df)
count_dict, strand_count_dict = get_min_counts(amin=0, amax=101, astep=1, dmin=200, dmax=1400, dstep=400, modelrefcount=refmean, bon=28119) ##CHANGE: variables and for bonf correction
plot_count_thresh_dict(count_dict, strand_count_dict, name='model_alt_nuc_simulation_minVAF_plot.png')    

