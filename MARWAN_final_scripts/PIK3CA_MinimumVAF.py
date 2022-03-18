import pandas as pd
import numpy as np
import pysam
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
pd.set_option('display.max_columns', None)
plt.style.use('ggplot')

# Read in the fasta and model file

fastafile = 'hg19.fa'
posmodelfile = 'OGPB_model.txt'

# Append a column to the dataframe with the reference allele
def get_refs(df, fasta):
    fa = pysam.FastaFile(fasta)
    for index, row in df.iterrows():
        chromo = index[0]
        pos = index[1]
        ref = fa.fetch(chromo, pos - 1, pos)
        df.loc[index, 'Ref'] = ref
    return df

df = pd.read_csv(posmodelfile, sep='\t')
df = df[df['Chromosome'] == 'chr3'] # This is positions 178,916,614 - 178,952,152 (Just PIK3CA [which is from 178,865,902-178,957,881])

# Read in the vcf file
vcffile = 'PIK3CA_varsifter_loc.txt.uniq'
vcf = pd.read_csv(vcffile, sep='\t', header=None).values

# Remove the intersected positions and append reference positions to the dataframe
df2 = df[df['Position'] != vcf[0,1]]
for i in range(len(vcf)):
    df2 = df2[df2['Position'] != vcf[i,1]]
df2.set_index(['Chromosome', 'Position'], inplace=True)
df2 = get_refs(df2, fastafile)

def get_min_VAF_for_positions(model, bon, depths):
    df = model
    nucs = ['A', 'C', 'G', 'T']
    count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    strand_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    # iterate through the 4 nucleotides A,C,G,T
    for nuc in nucs:
        #create an array called alts that contains the other 3 nucleotides. Ex: If on "A", then alts has ["C","G","T"]
        alts = [alt_allele for alt_allele in nucs if alt_allele != nuc]
        # subset the original dataframe for only positions with a refernce that is the same as the current nucleotide. Ex: If on "A", then only look at posiitons where "A" is the reference
        nuc_df = df[df['Ref'] == nuc]
        # Then iterate through the 3 possible alterate nucleotides
        for alt in alts:
            # For a single chosen alt, then iterate through the possible depths (in this case just 600X)
            for depth in depths:
                # Create an empty list
                strand_count_min_list = []
                # At the selected alternate allele and depth of coverage, iterate through the nuc_df dataframe (which is subsetted to contain only positions with the reference allele of interest)
                for index, row in nuc_df.iterrows():
                    # Initiate the alt counter at 1
                    alt_counter = 1
                    # Initiate threshold status to False - We haven't hit significance yet!
                    thresh = False
                    # While not at the threshold:
                    while not thresh:
                        ref_depth = int(round(depth / 2.))
                        for_count = int(round(alt_counter / 2.))
                        rev_count = (alt_counter - for_count)
                        altIDfor = alt+'_for'
                        altIDrev = alt+'_rev'
                        refIDfor = nuc + '_for'
                        refIDrev = nuc + '_rev'
                        for_table = np.array([[for_count, row[altIDfor]], [(ref_depth - for_count), row[refIDfor]]])
                        rev_table = np.array([[rev_count, row[altIDrev]], [(ref_depth - rev_count), row[refIDrev]]])
                        # Calculate forward strand p value
                        for_pval = stats.fisher_exact(for_table, alternative='greater')[1]
                        # Calculate reverse strand p value
                        rev_pval = stats.fisher_exact(rev_table, alternative='greater')[1]
                        # Combine the p values
                        strand_pval = stats.combine_pvalues([for_pval, rev_pval], method='stouffer')[1]
                        # Check if p value is significant (with the bonferoni correction)
                        if strand_pval < (0.05 / bon):
                            # If it is, then set the threshold value = True
                            thresh = True
                            # Set the vaf = number of alt alleles / depth
                            vaf = alt_counter/float(depth)
                            # Append that VAF value to the strand_count_min_list we created
                            strand_count_min_list.append(vaf)
                        # Otherwise, if the pval is not significant, then increase the alt counter and try again
                        else:
                            alt_counter += 1
                strand_count_dict[nuc][alt][depth] = strand_count_min_list
    return strand_count_dict

strand_min_dict = get_min_VAF_for_positions(model=df2, bon=28119, depths=[600])

# Figure 5a

C_data = [strand_min_dict['C']['A'][600],strand_min_dict['C']['G'][600],strand_min_dict['C']['T'][600]]
C_A = C_data[0]
C_G = C_data[1]
C_T = C_data[2]

x_data = np.arange(1, len(C_A) + 1, 1)

fig, ax = plt.subplots(3, sharey=True,figsize=(13.33, 7.5))

ax[0].plot(x_data,C_A, color='blue', linewidth=1, label='A')
ax[1].plot(x_data,C_G, color='green', linewidth=1, label='G')
ax[2].plot(x_data,C_T, color='goldenrod', linewidth=1, label='T')

for j in range(len(ax)):
    ax[j].tick_params(axis='x', colors='black')
    ax[j].tick_params(axis='y', colors='black')
    ax[j].set_ylim(0, 0.03)
    ax[j].set_xlim(0, len(x_data) + 1)
    ax[j].set_yticks([0.00,0.01,0.02,0.03])
    leg = ax[j].legend(prop={'family': 'Arial'}, loc=1)
    for line in leg.get_lines():
        line.set_linewidth(3.0)
ax[2].set_xlabel('Ordered reference allele position', color='black', fontname='Arial', fontsize=22, labelpad=15)
ax[1].set_ylabel('Minimum VAF', color='black', fontname='Arial', fontsize=22, labelpad=15)

plt.savefig('Figure_5A_test.png', dpi=600)

# Figure 5b for Cytosine

C_data = [strand_min_dict['C']['A'][600],strand_min_dict['C']['G'][600],strand_min_dict['C']['T'][600]]

fig = plt.figure(1, figsize=(6, 6))
ax = fig.add_subplot(111)
boxcolors = ['blue','green','goldenrod']
whiskercolors = ['blue','blue','green','green','goldenrod','goldenrod']

bp = ax.boxplot(C_data,patch_artist=True,widths=0.4)

for i,box in enumerate(bp['boxes']):
    # change outline color
    box.set( color=boxcolors[i], linewidth=2)
    # change fill color
    box.set( facecolor = 'None' )

## change color and linewidth of the whiskers
for j,whisker in enumerate(bp['whiskers']):
    whisker.set(color=whiskercolors[j], linewidth=2)

## change color and linewidth of the caps
for k,cap in enumerate(bp['caps']):
    cap.set(color=whiskercolors[k], linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='k', linewidth=2)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='k' ,linewidth=2)

ax.set_ylim(0)
ax.set_xticklabels(['A', 'G', 'T'])
ax.set_ylabel('Minimum VAF',color='black', fontname='Arial', fontsize=16)
ax.set_xlabel('Alternate nucleotide',color='black', fontname='Arial', fontsize=16)
ax.set_title('Cytosine Reference Allele',color='black', fontname='Arial', fontsize=16)
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

plt.savefig('Cytosine_Reference_Boxplot', dpi=600)

# For Figure 5b for Guanine

G_data = [strand_min_dict['G']['A'][600],strand_min_dict['G']['C'][600],strand_min_dict['G']['T'][600]]

fig = plt.figure(2, figsize=(6, 6))
ax = fig.add_subplot(111)
boxcolors = ['blue','red','goldenrod']
whiskercolors = ['blue','blue','red','red','goldenrod','goldenrod']

bp = ax.boxplot(G_data,patch_artist=True,widths=0.4)

for i,box in enumerate(bp['boxes']):
    # change outline color
    box.set( color=boxcolors[i], linewidth=2)
    # change fill color
    box.set( facecolor = 'None' )

## change color and linewidth of the whiskers
for j,whisker in enumerate(bp['whiskers']):
    whisker.set(color=whiskercolors[j], linewidth=2)

## change color and linewidth of the caps
for k,cap in enumerate(bp['caps']):
    cap.set(color=whiskercolors[k], linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='k', linewidth=2)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='k',linewidth=2)

ax.set_ylim(0)
ax.set_xticklabels(['A', 'C', 'T'])
ax.set_ylabel('Minimum VAF',color='black', fontname='Arial', fontsize=16)
ax.set_xlabel('Alternate nucleotide',color='black', fontname='Arial', fontsize=16)
ax.set_title('Guanine Reference Allele',color='black', fontname='Arial', fontsize=16)
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

plt.savefig('Guanine_Reference_Boxplot', dpi=600)

# For Figure 5b for Adenine

A_data = [strand_min_dict['A']['C'][600],strand_min_dict['A']['G'][600],strand_min_dict['A']['T'][600]]

fig = plt.figure(3, figsize=(6, 6))
ax = fig.add_subplot(111)
boxcolors = ['red','green','goldenrod']
whiskercolors = ['red','red','green','green','goldenrod','goldenrod']

bp = ax.boxplot(A_data,patch_artist=True,widths=0.4)

for i,box in enumerate(bp['boxes']):
    # change outline color
    box.set( color=boxcolors[i], linewidth=2)
    # change fill color
    box.set( facecolor = 'None' )

## change color and linewidth of the whiskers
for j,whisker in enumerate(bp['whiskers']):
    whisker.set(color=whiskercolors[j], linewidth=2)

## change color and linewidth of the caps
for k,cap in enumerate(bp['caps']):
    cap.set(color=whiskercolors[k], linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='k', linewidth=2)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='k' ,linewidth=2)

ax.set_ylim(0)
ax.set_xticklabels(['C', 'G', 'T'])
ax.set_ylabel('Minimum VAF',color='black', fontname='Arial', fontsize=16)
ax.set_xlabel('Alternate nucleotide',color='black', fontname='Arial', fontsize=16)
ax.set_title('Adenine Reference Allele',color='black', fontname='Arial', fontsize=16)
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

plt.savefig('Adenine_Reference_Boxplot', dpi=600)

# For Figure 5b for Thymine

T_data = [strand_min_dict['T']['A'][600],strand_min_dict['T']['C'][600],strand_min_dict['T']['G'][600]]

fig = plt.figure(4, figsize=(6, 6))
ax = fig.add_subplot(111)
boxcolors = ['blue','red','green']
whiskercolors = ['blue','blue','red','red','green','green']

bp = ax.boxplot(T_data,patch_artist=True,widths=0.4)

for i,box in enumerate(bp['boxes']):
    # change outline color
    box.set( color=boxcolors[i], linewidth=2)
    # change fill color
    box.set( facecolor = 'None' )

## change color and linewidth of the whiskers
for j,whisker in enumerate(bp['whiskers']):
    whisker.set(color=whiskercolors[j], linewidth=2)

## change color and linewidth of the caps
for k,cap in enumerate(bp['caps']):
    cap.set(color=whiskercolors[k], linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='k', linewidth=2)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='k' ,linewidth=2)

ax.set_ylim(0)
ax.set_xticklabels(['A', 'C', 'G'])
ax.set_ylabel('Minimum VAF',color='black', fontname='Arial', fontsize=16)
ax.set_xlabel('Alternate nucleotde',color='black', fontname='Arial', fontsize=16)
ax.set_title('Thymine Reference Allele',color='black', fontname='Arial', fontsize=16)
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

plt.savefig('Thymine_Reference_Boxplot', dpi=600)

# Data for Supplementary Table 4

nucs = ['A', 'C', 'G', 'T']
alts = ['A', 'C', 'G', 'T']
for nuc in nucs:
    for alt in alts:
        print(nuc,">",alt,np.mean(strand_min_dict[nuc][alt][600]))


