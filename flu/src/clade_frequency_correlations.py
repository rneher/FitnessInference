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
#
################################################################################

import matplotlib
matplotlib.use('pdf')
import sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
import test_flu_prediction as test_flu
import predict_flu as flu
import tree_utils
from Bio import Phylo,AlignIO,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date

params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': False}
plt.rcParams.update(params)
figure_folder = '../figures_ms/'
analysis_folder = '../analysis_may_feb/'

def clade_frequency_correlations_func(params):
    ssfactor = 1.0
    year=params.year
    prediction_regions = params.pred.split(',')
    test_regions = params.test.split(',')
    flutype = params.flutype
    aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'

    if flutype.startswith('H3N2'):
        cds = {'begin':0, 'end':987, 'pad':0}
    else:
        cds = {'begin':0, 'end':300, 'pad':0}

    base_name = flutype+'_'+str(year)+'_trees_pred_'+'_'.join(prediction_regions)+'_test_'+'_'.join(test_regions)
    base_name=base_name.replace('north america', 'NA')
    print base_name
    name_mod = 'ssize_'+str(params.sample_size)+'_b_'+str(params.boost)+'_d_'+str(params.dscale)+'_D_'+str(params.diffusion)

    # open annotations file
    with open('../data/'+flutype+'_annotations.pickle', 'r') as infile:
        annotation = pickle.load(infile)
    outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')

    # define prediction and test data sets by choosing start and stop dates 
    # begin previous year on may 1st, end this year on Feb 28
    prediction_set={'start': date(year-1, 5,1), 'stop':date(year, 2,28),
                    'regions':prediction_regions, 'sample_size':params.sample_size}
    # test data start Oct 1st, end following year end of March
    test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                'regions':test_regions, 'sample_size':params.sample_size}


    # define time intervals to estimate changing clade frequencies.
    bin_dt = 105
    tbins_pred = [ date.fromordinal(prediction_set['stop'].toordinal()-ii*bin_dt) for ii in range(
            (prediction_set['stop'].toordinal()-prediction_set['start'].toordinal())//bin_dt,-1,-1)]

    tbins_eval = [ prediction_set['stop'], test_set['stop']]

    # set up the prediction and pass all parameters to the wrapper function
    prediction = test_flu.predict(aln_fname, outgroup, annotation,\
                            ['mean_fitness', 'expansion_score', flu.combined_ranking_internal, flu.combined_ranking_external],
                            prediction_set, cds, time_bins = tbins_pred, subsample_factor = ssfactor, boost = params.boost,
                            eps_branch_length = params.eps_branch, dscale = params.dscale, D=params.diffusion, collapse = params.collapse)
    test_data = flu.flu_alignment(aln_fname,outgroup,annotation,cds=cds,
                             criteria = [[test_set['start'], test_set['stop'],
                            [reg],test_set['sample_size']] for reg in test_set['regions']], collapse=params.collapse)

    # combine the test data, the prediction data and possible other sequences
    # and build a tree
    combined_data = test_flu.make_combined_data(prediction, test_data, collapse = params.collapse)
    combined_tree = flu.flu_ranking(combined_data, time_bins = tbins_eval, pseudo_count = 0)
    combined_tree.expansion_score()

    tree_utils.find_internal_nodes(prediction.T,combined_tree.T)

    fitness_and_freqs = []
    for node in prediction.non_terminals:
        fitness_and_freqs.append([node.mean_fitness]+list(node.mirror_node.temporal_frequency))

    fitness_and_freqs = np.array(fitness_and_freqs)
    return fitness_and_freqs

class specs(object):
    def __init__(self, year = 1995, dscale = 5.0, eps_branch = 1e-5, flutype='H3N2', 
                 boost=0.0, collapse =False , diffusion=0.5, pred= 'asia,north america',
                 test = 'asia,north america', sample_size = 100):
        self.year=year
        self.dscale=dscale
        self.eps_branch=eps_branch
        self.flutype=flutype
        self.diffusion=diffusion
        self.sample_size=sample_size
        self.pred=pred
        self.test=test
        self.boost=boost
        self.collapse=collapse

def clade_frequency_correlations_all_years(dscale = 5.0, years = range(1995,2014)):
    fitness_and_freqs = []
    for year in years:
        params = specs(year=year, collapse=True, sample_size=200)
        fitness_and_freqs.append(clade_frequency_correlations_func(params))
    return fitness_and_freqs

def plot_clade_frequencies_correlations(fitness_and_freqs):
    all_freqs = []
    plt.figure()
    for faf in fitness_and_freqs:
        tmp = faf[(faf[:,1]<0.8)*(faf[:,1]>0.005)]
        rank_fit = stats.rankdata(tmp[:,0])
        rank_fit/=rank_fit.max()
        rank_freq = stats.rankdata(tmp[:,2]/tmp[:,1])
        rank_freq/=rank_freq.max()
        all_freqs.extend([(a,b) for a,b in zip(rank_fit,rank_freq)])
        plt.scatter(rank_fit, rank_freq)
        print stats.spearmanr(tmp[:,0],tmp[:,2]/tmp[:,1])

    plt.xlabel('fitness rank')
    plt.ylabel('growth rank')
    all_freqs =  np.array(all_freqs)
    print 'all', stats.spearmanr(all_freqs[:,0], all_freqs[:,1])
    return all_freqs


if __name__=="__main__":

    #parse the command line arguments
    parser = argparse.ArgumentParser(description="predict the flu")
    parser.add_argument('--flutype', type = str, default = 'H3N2', help = 'influenza strain')
    parser.add_argument('--year', default=-1, type=int, help='year/(year+1) season which is to be predicted')
    parser.add_argument('--sample_size', default=100, type=int, 
                        help='max number of sequences to use')
    parser.add_argument('--pred', default='asia,north america', type=str, help='regions to predict')
    parser.add_argument('--test', default='asia,north america', type=str, help='regions to test')
    parser.add_argument('--boost', default=0.0, type=float, help='fitness increment for cluster mutation')
    parser.add_argument('--eps_branch', default=1e-5, type=float, help='minimal branch length for inference')
    parser.add_argument('--diffusion', default=0.5, type=float, help='Fitness diffusion coefficient')
    parser.add_argument('--dscale', default=5.0, type=float, help='scale factor for time scale')
    parser.add_argument('--collapse', const = True, default=False, nargs='?', 
                        help='collapse internal branches with identical sequences')

    params=parser.parse_args()
    if params.year>0:
        faf = clade_frequency_correlations_func(params)
        plot_clade_frequencies_correlations([faf])
    else:
        faf = clade_frequency_correlations_all_years(dscale=params.dscale, years = range(2003,2014))
        normed_freqs = plot_clade_frequencies_correlations(faf)
        plt.savefig(figure_folder+'fitness_frequency_correlations_years_2003_2014_gamma_'+str(params.dscale)+'.pdf')
        with open(analysis_folder+'fitness_frequency_correlations_years_2003_2014_gamma_'+str(params.dscale)+'.pickle', 'w') as outfile:
            pickle.dump(normed_freqs, outfile)
        
        
