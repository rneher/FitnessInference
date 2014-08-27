#!/ebio/ag-neher/share/programs/bin/python2.7
################################################################################
#
# author: Richard Neher
# email:  richard.neher@tuebingen.mpg.de
#
# Reference: Richard A. Neher, Colin A Russell, Boris I Shraiman. 
#            "Predicting evolution from the shape of genealogical trees"
#
################################################################################
#
#script that analyzes  repeatedly a set of aligned sequences of the HA1 domain of influenza
#and attempts to predicts a strain close the population of the following year.
#We operate on historical data, hence the prediction is tested right away on
#next years strain collection. It does not plot, but saves predictions on repeated 
#subsamples to file.
#It takes arguments and alignment, a year, and sample size as optional arguments.
#more detailed parameters can only be set by hand as of now.
#
################################################################################

import sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
sys.path.append('../../src')
import test_flu_prediction as test_flu
import predict_flu as flu
import tree_utils
from Bio import Phylo,AlignIO,SeqIO
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date

analysis_folder = '../analysis_may_feb/'

#parse the command line arguments
parser = argparse.ArgumentParser(description="predict the flu")
parser.add_argument('--flutype', type = str, default = 'H3N2', help = 'influenza strain')
parser.add_argument('--year', default=2000, type=int, help='year/(year+1) season which is to be predicted')
parser.add_argument('--sample_size', default=100, type=int, 
                    help='max number of sequences to use')
parser.add_argument('--pred', default='asia,north america', type=str, help='regions to predict')
parser.add_argument('--test', default='asia,north america', type=str, help='regions to test')
parser.add_argument('--boost', default=0.0, type=float, help='fitness increment for cluster mutation')
parser.add_argument('--eps_branch', default=1e-5, type=float, help='minimal branch length for inference')
parser.add_argument('--diffusion', default=0.5, type=float, help='Fitness diffusion coefficient')
parser.add_argument('--dscale', default=5.0, type=float, help='scale factor for time scale')
parser.add_argument('--collapse', const = True, default=False, nargs='?', help='collapse internal branches with identical sequences')
parser.add_argument('--nreps', default=10, type=int, help='number of repetitions')

ssfactor = 0.7
params=parser.parse_args()
flutype = params.flutype
year=params.year
prediction_regions = params.pred.split(',')
test_regions = params.test.split(',')
pseudo_count = 5

aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'

if flutype.startswith('H3N2'):
    cds = {'begin':0, 'end':987, 'pad':0}
else:
    cds = {'begin':0, 'end':300, 'pad':0}
    
base_name = flutype+'_'+str(year)+'_trees_asia_NA_'
name_mod = 'ssize_'+str(params.sample_size)+'_b_'+str(params.boost)+'_d_'\
    +str(params.dscale)+'_D_'+str(params.diffusion)
if params.collapse: name_mod+='_collapse'

# open annotations file
with open('../data/'+flutype+'_annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)
outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')

# define prediction and test data sets by choosing start and stop dates 
if "oceania" in test_regions:
    prediction_set={'start': date(year-2, 10,1), 'stop':date(year-1, 9,30),
                    'regions':prediction_regions, 'sample_size':params.sample_size}
    test_set = {'start':date(year, 3,1), 'stop':date(year, 9,30),
                'regions':test_regions, 'sample_size':params.sample_size}
else:
    # begin previous year on may 1st, end this year on Feb 28
    prediction_set={'start': date(year-1,5,1), 'stop':date(year, 2,28),
                    'regions':prediction_regions, 'sample_size':params.sample_size}
    test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                'regions':test_regions, 'sample_size':params.sample_size}

# define 105 day intervals to estimate changing clade frequencies.
# chosen to have 3 intervals between May and Feb
bin_dt = 105
tbins = [ date.fromordinal(prediction_set['stop'].toordinal()-ii*bin_dt) for ii in range(
        (prediction_set['stop'].toordinal()-prediction_set['start'].toordinal())//bin_dt,-1,-1)]
    
# allocate arrays to save the predictions
nuc_dist_array = np.zeros((params.nreps, 10))
epi_dist_array = np.zeros((params.nreps, 10))
top_strains = []
for ii in xrange(params.nreps):
    # set up the prediction and pass all parameters to the wrapper function
    prediction = test_flu.predict(aln_fname, outgroup, annotation,\
                            ['mean_fitness', 'expansion_score', 'depth', flu.combined_ranking_internal, flu.combined_ranking_external],
                            prediction_set, cds, time_bins = tbins, subsample_factor = ssfactor, boost = params.boost,
                            eps_branch_length = params.eps_branch, collapse = params.collapse, dscale = params.dscale, D=params.diffusion, pseudo_count = pseudo_count)

    # define the methodes for which the predictions are to be evaluated
    methods = [('mean_fitness', '_ext', prediction.terminals),
                ('mean_fitness', '_int', prediction.non_terminals),
                ('expansion_score', '_int', prediction.non_terminals),
                ('expansion_fitness', '', prediction.non_terminals),
                ('time_fitness', '', prediction.terminals),
                ('ladder_rank', '', prediction.terminals),
                ('date', '', prediction.terminals)]
    distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, test_set, aln_fname, outgroup, annotation,cds=cds, flutype=flutype)
    nuc_dist_array[ii,:] = [distances['average'],distances['minimal'],distances['L&L']]\
                    +[distances[m[0]+m[1]] for m in methods]
    epi_dist_array[ii,:] = [distances_epi['average'],distances_epi['minimal'],distances_epi['L&L']]\
                    +[distances_epi[m[0]+m[1]] for m in methods]
    # memorize the strain predicted best
    top_strains.append(prediction.best_node(method = 'mean_fitness', nodes = prediction.terminals))

# save the results to file
fname_nuc = '_'.join(map(str,[analysis_folder+flutype,'year',year, 'pred']+ 
                     ["-".join(prediction_regions)]+['test']+
                     ['-'.join(test_regions),name_mod])).replace(' ', '')+'_nuc.dat'
#if file does not exist, create and write header
if not os.path.isfile(fname_nuc):
    with open(fname_nuc, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_nuc, 'a') as outfile:
    np.savetxt(outfile, nuc_dist_array)


# save the results to file
fname_nuc = '_'.join(map(str,[analysis_folder+flutype,'year',year, 'pred']+ 
                     ["-".join(prediction_regions)]+['test']+
                     ['-'.join(test_regions),name_mod])).replace(' ', '')+'_epi.dat'
#if file does not exist, create and write header
if not os.path.isfile(fname_nuc):
    with open(fname_nuc, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_nuc, 'a') as outfile:
    np.savetxt(outfile, epi_dist_array)


# write the best nodes to file:
fname_strains = '_'.join(map(str,[analysis_folder+flutype,'year',year, 'topstrains']+ 
                     ["-".join(prediction_regions)]+['test']+
                     ['-'.join(test_regions),name_mod])).replace(' ', '')+'.dat'
with open(fname_strains, 'a') as outfile:
    for strain in top_strains:
        outfile.write(strain.name+'\t'+str(strain.date)+'\n')

