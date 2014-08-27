'''
This script analyze previously computed tables of the fitness differences along branches,
the number of synonymous, non-synonymous, epitope mutations, adn koel mutations
These files are computed by analyze_HA1.py
'''

import numpy as np
import glob, pickle, argparse
from scipy import stats

parser = argparse.ArgumentParser(description="analyze mutation prevalence on branches")
parser.add_argument('--flutype', type = str, default = 'H3N2', help = 'influenza strain')
parser.add_argument('--sample_size', default=100, type=int, 
                    help='max number of sequences to use')
parser.add_argument('--pred', default='asia,north america', type=str, help='regions to predict')
parser.add_argument('--test', default='asia,north america', type=str, help='regions to test')
parser.add_argument('--boost', default=0.0, type=float, help='fitness increment for cluster mutation')
parser.add_argument('--diffusion', default=0.5, type=float, help='Fitness diffusion coefficient')
parser.add_argument('--dscale', default=5.0, type=float, help='scale factor for time scale')
parser.add_argument('--collapse', const = True, default=False, nargs='?', help='collapse internal branches with identical sequences')

params=parser.parse_args()
prediction_regions = params.pred.split(',')
test_regions = params.test.split(',')

# construct file mask to read precomputed files
analysis_folder = '../analysis_may_feb/'
figure_folder = '../figures_ms/'
base_name = params.flutype+'_????_trees_pred_'+'_'.join(prediction_regions)\
            +'_test_'+'_'.join(test_regions)
base_name=base_name.replace('north america', 'NA')
print base_name
name_mod = 'ssize_'+str(params.sample_size)+'_b_'+str(params.boost)+'_d_'+str(params.dscale)\
            +'_D_'+str(params.diffusion)
if params.collapse: name_mod+='_collapse'

#make a list of all relevant files
flist = glob.glob(analysis_folder+base_name+'*_dfit_mutations_*'+name_mod+'.dat')

# load and make a long list of all branches
dfit_all = []
for fname in flist:
    dfit_all.extend(np.loadtxt(fname))

# case to array to enable slicing
dfit_all=np.array(dfit_all)
# the second column contains the terminal/non-terminal designation. 1==terminal
# select only non-terminal branches
non_terminals = dfit_all[:,1]==0
dfit_non_terminals = dfit_all[non_terminals,:]

# loop over different slices of the distribution of delta fitness
# and summarize mutation counts of different categories
with open(figure_folder+params.flutype+'_dfit_mutations_'+name_mod+'.dat', 'w') as outfile:
    prev_threshold = -10000
    total_mut_numbers = []
    delim = ' & '
    outstr = delim.join(map(str,['perc', '\# non-syn','\# syn', '\# epi', '\# Koel', 'dn/ds', 'epi/ds', 'Koel/ds']))
    print outstr
    outfile.write(outstr+'\n')
    for perc_cutoff in [25,50,75,100]:
        threshold = stats.scoreatpercentile(dfit_non_terminals[:,0], perc_cutoff)
        # pull out indices belonging to the following delta fitness window
        # prev_threshold <= dfit < threshold
        selected_fit = (dfit_non_terminals[:,0]>=prev_threshold) \
                    * (dfit_non_terminals[:,0]< threshold)
    
        # look at only these indices
        dn = np.sum(dfit_non_terminals[selected_fit,3])
        ds = np.sum(dfit_non_terminals[selected_fit,2]) - dn
        epi = np.sum(dfit_non_terminals[selected_fit,4])
        cluster = np.sum(dfit_non_terminals[selected_fit,5])
        total_mut_numbers.append([dn,ds,epi,cluster])
    
        # print resulting numbers
        outstr = delim.join(map(str,[perc_cutoff, dn,ds, epi, cluster, round(dn/ds,3), round(epi/ds,3), round(cluster/ds,3)]))
        print outstr
        outfile.write(outstr+'\n')
    
        # update threshold
        prev_threshold = threshold
    
    # cast to array to enable slicing
    total_mut_numbers=np.array(total_mut_numbers)
    # make a table with Fisher-Exact test pvals
    outstr = '\nupper vs lower delta fitness quartile:'
    print outstr
    outfile.write(outstr+'\n')
    for ci, cname in [(0,'dn vs ds'), (2,'epi vs ds'), (3,'Koel vs ds')]: 
        fs = stats.fisher_exact([[total_mut_numbers[0,1], total_mut_numbers[0,ci]],
                                [total_mut_numbers[-1,1], total_mut_numbers[-1,ci]]])
        outstr = delim.join(map(str,[cname, 'odds:',round(fs[0],3), 'pval:',round(fs[1],6)]))
        print outstr
        outfile.write(outstr+'\n')
    
    outstr = '\nupper vs lower delta fitness quartile:'
    print outstr
    outfile.write(outstr+'\n')
    for ci, cname in [(2,'epi vs dn'), (3,'Koel vs dn')]: 
        fs = stats.fisher_exact([[total_mut_numbers[0,0], total_mut_numbers[0,ci]],
                                [total_mut_numbers[-1,0], total_mut_numbers[-1,ci]]])
        outstr = delim.join(map(str,[cname, 'odds:',round(fs[0],3), 'pval:',round(fs[1],6)]))
        print outstr
        outfile.write(outstr+'\n')

