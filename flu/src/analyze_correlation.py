import sys
sys.path.append('../../prediction_src')
import test_flu_prediction as flu
from Bio import Phylo,AlignIO,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr

corrfunc = spearmanr
#corrfunc = pearsonr

ssfactor=1.0
prediction_regions = ('asia', 'north america')
test_regions = prediction_regions
sample_size = 200
boost = 0.0
eps_branch_length = 1e-5
dscale = 5.0
collapse = False
flutype = 'H3N2'
D = 0.5
pseudo_count = 5
min_frac = 0.05
max_frac = 0.5

aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'

if flutype.startswith('H3N2'):
    cds = {'begin':0, 'end':987, 'pad':0}
else:
    cds = {'begin':0, 'end':300, 'pad':0}
    
if os.path.isfile('../data/'+flutype+'_L_L_predictions.pickle'):
    with open('../data/'+flutype+'_L_L_predictions.pickle') as infile:
        laessig_prediction = pickle.load(infile)

# open annotations file
with open('../data/'+flutype+'_annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)
outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')

bin_dt = 105  #time bins in days. 3*105 = 315 days approx 10 month
years = range(1995,2012)
predictions = {}
for year in years:
    if "oceania" in test_regions:
        prediction_set={'start': date(year-2, 10,1), 'stop':date(year-1, 9,30),
                        'regions':prediction_regions, 'sample_size':sample_size}
        test_set = {'start':date(year, 3,1), 'stop':date(year, 9,30),
                    'regions':test_regions, 'sample_size':sample_size}
    else:
        prediction_set={'start': date(year-1, 5,1), 'stop':date(year, 2,28),
                        'regions':prediction_regions, 'sample_size':sample_size}
        test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                    'regions':test_regions, 'sample_size':sample_size}

        tbins = [ date.fromordinal(prediction_set['stop'].toordinal()-ii*bin_dt) for ii in range(
        (prediction_set['stop'].toordinal()-prediction_set['start'].toordinal())//bin_dt,-1,-1)]

    
    predictions[year] = flu.predict(aln_fname, outgroup, annotation,\
                            ['mean_fitness', 'expansion_score', flu.combined_ranking_internal, flu.combined_ranking_external],
                            prediction_set, cds, time_bins = tbins, subsample_factor = ssfactor, boost = boost,
                            eps_branch_length = 1e-5, collapse = collapse, dscale = dscale, D=D, pseudo_count = pseudo_count)

method_correlation = defaultdict(list)
for year in years:
    min_size = len(predictions[year].data.aln)*min_frac
    max_size = len(predictions[year].data.aln)*max_frac

    terminal_methods = np.array([ (a.mean_fitness, a.date.toordinal(), a.ladder_rank, -a.branch_length) for a in predictions[year].T.get_terminals()]).T
    method_correlation['ext_fit_date'].append(spearmanr(terminal_methods[0], terminal_methods[1])+ pearsonr(terminal_methods[0], terminal_methods[1]))
    method_correlation['ext_fit_ladder'].append(spearmanr(terminal_methods[0], terminal_methods[2])+pearsonr(terminal_methods[0], terminal_methods[2]))
    method_correlation['ext_fit_negbranchlength'].append(spearmanr(terminal_methods[0], terminal_methods[3])+ pearsonr(terminal_methods[0], terminal_methods[3]))

    nonterminal_methods = np.array([ (a.mean_fitness, a.expansion_score, a.count_terminals()) 
                                     for a in predictions[year].T.get_nonterminals() 
                                     if a.count_terminals()>min_size and a.count_terminals()<max_size]).T
    method_correlation['int_fit_expansion'].append(spearmanr(nonterminal_methods[0], nonterminal_methods[1])+ pearsonr(nonterminal_methods[0], nonterminal_methods[1]))
    method_correlation['int_fit_count'].append(spearmanr(nonterminal_methods[0], nonterminal_methods[2])+pearsonr(nonterminal_methods[0], nonterminal_methods[2]))


for m in method_correlation:
    method_correlation[m]=np.array(method_correlation[m])
    plt.plot(years,method_correlation[m][:,0], label=m)

plt.legend(loc='lower center')
plt.xlabel('years')
plt.ylabel('rank correlation with fitness estimate')
plt.savefig('../figures/correlation_of_methods_'+'_'.join(map(str, ['dscale', dscale, 'D', D]))+'.pdf')
