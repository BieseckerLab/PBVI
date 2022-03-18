import os
import pandas as pd
import numpy as np
from scipy import stats
import glob
from collections import defaultdict
from functools import partial
import pysam
import matplotlib.pyplot as plt
plt.style.use('ggplot')
pd.set_option('display.max_columns', None)
import time

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

    #Frequency values are calculated and the alternate allele is determined
    #The following columns are added to the df: freqA	freqC	freqG	freqT	freqAltA	freqAltC	freqAltG	freqAltT	AltACount	AltCCount	AltGCount	AltTCount	freqAlt	Alt	AltCount

    return df

def get_alt_size_pvals(df, posmodel, size, CountThresh=0, alternative='two-sided', method='stouffer'):

    #Read in the df for the model
    posmodel_df = pd.read_csv(posmodel, sep='\t', index_col=[0, 1])
    df.set_index(['Chromosome', 'Position'], inplace=True)

    #FILTER#
    #Look at only positions with an alt allele count > threshold (default is 0 -- so anything with at least 1 alt allele)
    df = df[df['AltCount'] > CountThresh]

    #Join the pos_model df to df with statistics from the get_alt_stats function (adding the model info to the data frame)
    df = df.join(posmodel_df, lsuffix='_var', rsuffix='_posmodel')

    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_choices = [df['A_var'], df['C_var'], df['G_var'], df['T_var']]
    model_choices = [df['A_posmodel'], df['C_posmodel'], df['G_posmodel'], df['T_posmodel']]
    df['TestPosRef'] = np.select(ref_conditions, var_choices)
    df['TestPosAlt'] = np.select(alt_conditions, var_choices)
    df['ModelPosRef'] = np.select(ref_conditions, model_choices) * size
    df['ModelPosAlt'] = np.select(alt_conditions, model_choices) * size
    df.dropna(subset=['ModelPosRef'], inplace=True)
    df['ModelPosRef'] = df['ModelPosRef'].round(0).astype(int)
    df['ModelPosAlt'] = df['ModelPosAlt'].round(0).astype(int)
    # Calculate p value
    df['PosPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosAlt, r.ModelPosAlt], [r.TestPosRef, r.ModelPosRef]], alternative=alternative), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosPval'] = df['PosPval'].astype(float)

    # At the point, the following 5 columns have been added to the df: TestPosRef,TestPosAlt,ModelPosRef,ModelPosAlt,PosPval

    # The following block of code just splits up the count data for the "test" and "model" into forward and reverse strand specific

    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_for_choices = [df['A_for_var'], df['C_for_var'], df['G_for_var'], df['T_for_var']]
    model_for_choices = [df['A_for_posmodel'], df['C_for_posmodel'], df['G_for_posmodel'], df['T_for_posmodel']]
    df['TestPosForRef'] = np.select(ref_conditions, var_for_choices)
    df['TestPosForAlt'] = np.select(alt_conditions, var_for_choices)
    df['ModelPosForRef'] = np.select(ref_conditions, model_for_choices) * size
    df['ModelPosForAlt'] = np.select(alt_conditions, model_for_choices) * size
    df['ModelPosForRef'] = df['ModelPosForRef'].round(0).astype(int)
    df['ModelPosForAlt'] = df['ModelPosForAlt'].round(0).astype(int)
    # Calculate p value forward
    df['PosForPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosForAlt, r.ModelPosForAlt], [r.TestPosForRef, r.ModelPosForRef]], alternative=alternative), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosForPval'] = df['PosForPval'].astype(float)

    alt_conditions = [df['Alt'] == 'A', df['Alt'] == 'C', df['Alt'] == 'G', df['Alt'] == 'T']
    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    var_rev_choices = [df['A_rev_var'], df['C_rev_var'], df['G_rev_var'], df['T_rev_var']]
    model_rev_choices = [df['A_rev_posmodel'], df['C_rev_posmodel'], df['G_rev_posmodel'], df['T_rev_posmodel']]
    df['TestPosRevRef'] = np.select(ref_conditions, var_rev_choices)
    df['TestPosRevAlt'] = np.select(alt_conditions, var_rev_choices)
    df['ModelPosRevRef'] = np.select(ref_conditions, model_rev_choices) * size
    df['ModelPosRevAlt'] = np.select(alt_conditions, model_rev_choices) * size
    df['ModelPosRevRef'] = df['ModelPosRevRef'].round(0).astype(int)
    df['ModelPosRevAlt'] = df['ModelPosRevAlt'].round(0).astype(int)
    # Calculate p value reverse
    df['PosRevPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosRevAlt, r.ModelPosRevAlt], [r.TestPosRevRef, r.ModelPosRevRef]], alternative=alternative), axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['PosRevPval'] = df['PosRevPval'].astype(float)


    # At this point the following columns have been added to the df: TestPosForRef	TestPosForAlt	ModelPosForRef	ModelPosForAlt	PosForPval	TestPosRevRef	TestPosRevAlt	ModelPosRevRef	ModelPosRevAlt	PosRevPval

    # Calculate p value combined
    df['PosCombinePval'] = df.apply(lambda r: stats.combine_pvalues([r.PosForPval, r.PosRevPval], method=method), axis=1).str[1]

    # Added 1 more column: PosCombinePval

    ref_conditions = [df['Ref'] == 'A', df['Ref'] == 'C', df['Ref'] == 'G', df['Ref'] == 'T']
    ref_for_choices = [(df['A_for_var']), (df['C_for_var']), (df['G_for_var']), (df['T_for_var'])]
    df['TestPosForRef'] = np.select(ref_conditions, ref_for_choices)
    ref_rev_choices = [(df['A_rev_var']), (df['C_rev_var']), (df['G_rev_var']), (df['T_rev_var'])]
    df['TestPosRevRef'] = np.select(ref_conditions, ref_rev_choices)
    # Calculate p value (for strand bias)
    df['SBPval'] = df.apply(lambda r: stats.fisher_exact([[r.TestPosForAlt, r.TestPosRevAlt], [r.TestPosForRef, r.TestPosRevRef]]),axis=1).astype(str).str.split(',', 1, expand=True)[1].str[:-1]
    df['SBPval'] = df['SBPval'].astype(float)

    # Added final column: SBPval

    return df

def get_tp_fp_fn(df, model, thresh, n, truevars):

    sig_df = df[(df[model] * n) < thresh]
    
    #Count number of true positives, false positives, and false negatives
    total_tp = len(sig_df[sig_df.target == 1])
    total_fp = len(sig_df[sig_df.target == 0])
    total_fn = truevars - total_tp
    
    return total_tp, total_fp, total_fn

def get_model_f1(counts, model,varfiles_fn, bonf_correction):


    tps = []
    fps = []
    fns = []
    totalvars = []
    dfs = []
    for c,count in enumerate(counts):
        count_df = count.copy()
        count_df.reset_index(inplace=True, drop=False)

        ID = count_df.ID.iloc[0]
        print(ID)
        #This is the simulation output file of positions with simulated variants
        varfile = varfiles_fn[c]

        #Read in the varfile
        var_df = pd.read_csv(varfile, sep='\s+')
        #Only care about chr, pos, and alt allele
        var_df = var_df[['chromosome', 'position.1', 'alt_allele']]
        #Converted the var_df dataframe to a list of tuples
        target_list = list(var_df[['chromosome', 'position.1', 'alt_allele']].itertuples(index=False, name=None))

        #Add the target_list list of tuples as a column to the dataframe.
        count_df['var_tuple'] = list(count_df[['Chromosome', 'Position', 'Alt']].itertuples(index=False, name=None))

        #FILTER#
        #Only want to look at positions that have SBPval (2x2 fishers exact b/w TestPosForAlt,TestPosRevAlt ; TestPosForRef,TestPosRevRef) above bonferoni corrected. Want to only look at positions where there's no strand bias in the simulated BAM
        count_df = count_df[count_df['SBPval'] >= (0.05 / bonf_correction)]

        #Add target column. Assign a "1" if the position was in our original simulated positions list (from the BED-like varfile), and otherwise assign a 0.
        count_df['target'] = np.where(count_df['var_tuple'].isin(target_list), int(1), int(0))
        count_df.to_csv('marked.model.sigs.%s.txt' % ID, sep='\t')

        #Total number of positions in the BED-like variant file
        totaltargets = len(var_df)
        #Append this to a list. One element for each variant file, because each BAM file has a different set of simulated variants.
        totalvars.append(totaltargets)
        #Use the get_tp_fp_fn function to calculate true positives, false positives, and true negatives for the given BAM file of the loop
        tp, fp, fn = get_tp_fp_fn(df=count_df, model=model, thresh=0.05, n=bonf_correction, truevars=totaltargets)

        #Append the final count dataframe (i rows x 66 columns), true positives, false positives, and false negatives to 4 different lists
        dfs.append(count_df)
        tps.append(tp)
        fps.append(fp)
        fns.append(fn)


    #Sum the tps,fps,and fns across all tested BAM files
    tps = sum(tps)
    fps = sum(fps)
    fns = sum(fns)

    #Sum all the positions tested (Do I need to be concerned about duplicate positions?)
    gt = sum(totalvars)

    #Calculate precision
    if tps == 0 or fps == 0:
        precision = 0
    else:
        precision = float(tps) / (tps + fps)

    #Calculate recall (aka Sensitivity)
    if tps == 0 or fns == 0:
        recall = 0
    else:
        recall = float(tps) / (tps + fns)

    #Calculate F1 metric
    if precision == 0 or recall == 0:
        f1 = 0
    else:
        f1 = 1 / (((1 / precision) + (1 / recall)) / 2)

    return precision, recall, f1

def get_f1_results(bams, countfiles_fn, IDs_fn, pmodel, models, alternative='greater', minsize=0.1, maxsize=1.1, step=0.1):


    #Create a file with columns "size", "PosPval", "PosCombinePval" - for F1
    sizef = open('modelsize_f1_results_.'+str(minsize)+'_'+str(maxsize)+'_'+str(step)+'_'+'.txt', 'w')
    sizef.write('size\tPosPval\tPosCombinePval\n')
    size_dict = defaultdict(lambda: defaultdict(float))

    #Create an array where the elements are the proportion of the model size which will be tested
    model_sizes = np.arange(minsize, maxsize, step)

    #Create another file with columns "size", "PosPval", "PosCombinePval" - for Precision
    precf = open('modelsize_prec_results_.'+str(minsize)+'_'+str(maxsize)+'_'+str(step)+'_'+'.txt', 'w')
    precf.write('size\tPosPval\tPosCombinePval\n')
    prec_dict = defaultdict(lambda: defaultdict(float))

    #Create another file with columns "size", "PosPval", "PosCombinePval" - for Sensitivity
    recf = open('modelsize_rec_results_.'+str(minsize)+'_'+str(maxsize)+'_'+str(step)+'_'+'.txt', 'w')
    recf.write('size\tPosPval\tPosCombinePval\n')
    rec_dict = defaultdict(lambda: defaultdict(float))

    #Iterate through the model sizes
    for size in model_sizes:
        #Print which size we are currently testing
        print("Started analysis for size: " + str(size))
        #For each size, we will create a list called size_counts
        size_counts = []
        #Within each model size, we will iterate through all the bam files
        for q,bam in enumerate(bams):
            print("    Started analysis for bam file: " + bams[q])
            #Define the count file, which is created from the PBVI script
            countfile = countfiles_fn[q]
            #Read in the count file as a pandas dataframe (df)
            df = pd.read_csv(countfile, sep='\t')
            original_length = len(df)
            
            #Use the get_alt_stats to generate statistics for the initial dataframe and append those statistics to the dataframe
            alt_df = get_alt_stats(df)

            var_df = get_alt_size_pvals(alt_df, pmodel, size, 0, alternative, 'stouffer')

            #Get the ID for the bam file
            ID = IDs_fn[q]
            
            #Add the ID for the bam file to the ID column in the dataframe
            var_df['ID'] = ID

            #Append the newly created dataframe to the size_counts list. Each element in the list will be a dataframe for each bamfile
            size_counts.append(var_df)

        #After iterating through all the bam files, iterate through all the models we selected (just 2)
        #This portion of the code is calculating precision, recall, and f1 values for PosPval or PosCombinePval
        for model in models:
            #Print which model we are currently testing
            print("    Started analysis for: " + model)

            #Use get_model_f1 to calculate precision, recall, and f1 values for the given model and size
            prec, rec, f1 = get_model_f1(size_counts, model,varfiles,original_length)

            #Add these values to their respective dictionaries
            size_dict[size][model] = f1
            prec_dict[size][model]= prec
            rec_dict[size][model] = rec


        #After iterating through both models, then write dictionaries to the files we created at the very beginning. Repeat this for all sizes
        sizef.write(str(size)+'\t')
        sizef.write('\t'.join([str(size_dict[size][model]) for model in models]) + '\n')
        precf.write(str(size)+'\t')
        precf.write('\t'.join([str(prec_dict[size][model]) for model in models]) + '\n')
        recf.write(str(size)+'\t')
        recf.write('\t'.join([str(rec_dict[size][model]) for model in models]) + '\n')

    print("Finished analysis")

#First BAM

bamfiles = ['SAMPLE1.600X_simulation.bam']

varfiles = ['SAMPLE1.600X_simulation_output.txt']

countfiles = ['SAMPLE1.600X_PBVI_output_unannotated.txt']

IDs = ['SAMPLE1']

posmodelfile = 'OGPB_model.txt'

models = ['PosPval', 'PosCombinePval']

get_f1_results(bams=bamfiles, countfiles_fn=countfiles, IDs_fn=IDs, pmodel=posmodelfile,
               models=models, alternative='greater', minsize=0, maxsize=1.05, step=0.05)

sens_SAMPLE1 = pd.read_csv('modelsize_rec_results_.0_1.05_0.05_.txt', sep='\t')


#Second BAM

bamfiles = ['SAMPLE2.600X_simulation.bam']

varfiles = ['SAMPLE2.600X_simulation_output.txt']

countfiles = ['SAMPLE2.600X_PBVI_output_unannotated.txt']

IDs = ['SAMPLE2']

posmodelfile = 'OGPB_model.txt'

models = ['PosPval', 'PosCombinePval']

get_f1_results(bams=bamfiles, countfiles_fn=countfiles, IDs_fn=IDs, pmodel=posmodelfile,
               models=models, alternative='greater', minsize=0, maxsize=1.05, step=0.05)

sens_SAMPLE2 = pd.read_csv('modelsize_rec_results_.0_1.05_0.05_.txt', sep='\t')


#Third BAM

bamfiles = ['SAMPLE3.600X_simulation.bam']

varfiles = ['SAMPLE3.600X_simulation_output.txt']

countfiles = ['SAMPLE3.600X_PBVI_output_unannotated.txt']

IDs = ['SAMPLE3']

posmodelfile = 'OGPB_model.txt'

models = ['PosPval', 'PosCombinePval']

get_f1_results(bams=bamfiles, countfiles_fn=countfiles, IDs_fn=IDs, pmodel=posmodelfile,
               models=models, alternative='greater', minsize=0, maxsize=1.05, step=0.05)

sens_SAMPLE3 = pd.read_csv('modelsize_rec_results_.0_1.05_0.05_.txt', sep='\t')

# Average

sens_TOTAL = np.concatenate((sens_SAMPLE1.values[:,2].reshape(-1,1),sens_SAMPLE2.values[:,2].reshape(-1,1),sens_SAMPLE3.values[:,2].reshape(-1,1)),1)
sens_AVG = np.mean(sens_TOTAL,axis=1).reshape(-1,1)

# Figure 4a

sizesthing = np.linspace(0,1,len(sens_AVG)).reshape(-1,1)
plt.figure(figsize=(8, 7))
ax = plt.subplot(111)

plt.plot(sizesthing, sens_AVG, linewidth=3)
plt.ylim((0,1))
ax.set_ylabel('Sensitivity', color='black', fontname='Arial', fontsize=16)
ax.set_xlabel('Proportion of full model', color='black', fontname='Arial', fontsize=16)
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

plt.savefig('Figure_4A.png', dpi=600)
plt.show()

