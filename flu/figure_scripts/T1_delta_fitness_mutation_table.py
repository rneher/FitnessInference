'''
This script analyze previously computed tables of the fitness differences along branches,
the number of synonymous, non-synonymous, epitope mutations, adn koel mutations
These files are computed by analyze_HA1.py
'''
import sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
import test_flu_prediction as test_flu
import numpy as np
import glob, pickle
from scipy import stats

analysis_folder = test_flu.flu_analysis_folder
# parse the commandline arguments
parser = test_flu.make_flu_parser()
parser.add_argument('--tau', default = 1.0, type = float, help= 'memory time scale of the tree polarizer')
params=parser.parse_args()
params.year='????'

# get name snippets to link data files to requested parameters
base_name, name_mod = test_flu.get_fname(params)

figure_folder = '../figures_ms/'

#make a list of all relevant files
#fmask = analysis_folder+base_name+'_'+name_mod+'_dfit_mutations.dat'
fmask = analysis_folder+base_name+'_tau_'+str(params.tau)+'_polarizer_dfit_mutations.dat'
flist = glob.glob(fmask)

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
#with open(figure_folder+params.flutype+'_dfit_mutations_'+name_mod+'.dat', 'w') as outfile:
with open(figure_folder+params.flutype+'_dfit_mutations_tau_'+str(params.tau)+'_ssize_'+str(params.sample_size)+'.dat', 'w') as outfile:
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
    ORs = []
    for ci, cname in [(0,'dn vs ds'), (2,'epi vs ds'), (3,'Koel vs ds')]: 
        fs = stats.fisher_exact([[total_mut_numbers[0,1], total_mut_numbers[0,ci]],
                                [total_mut_numbers[-1,1], total_mut_numbers[-1,ci]]])
        ORs.append(fs)
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
        ORs.append(fs)
        print outstr
        outfile.write(outstr+'\n')

