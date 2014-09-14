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

file_formats = ['.pdf', '.svg']

plt.rcParams.update(test_flu.mpl_params)

figure_folder = '../figures_ms/'
analysis_folder = test_flu.flu_analysis_folder
mem_scale = 2.0**np.arange(-6,3)

def clade_frequency_correlations_func(params):
    # set up the prediction and pass all parameters to the wrapper function
    prediction = test_flu.predict_params(['mean_fitness'],
                                        params)


    # define the methodes for which the predictions are to be evaluated
    methods = [('mean_fitness', '_ext', prediction.terminals),
                ('mean_fitness', '_int', prediction.non_terminals)]
    distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)


    tbins_eval = [date(year=params.year, month = 6, day=1), date(year=params.year+1, month = 6, day=1)]
    combined_data = test_flu.make_combined_data(prediction, test_data, collapse = params.collapse)
    combined_tree = flu.flu_ranking(combined_data, time_bins = tbins_eval, pseudo_count = 0)
    combined_tree.expansion_score()

    tree_utils.find_internal_nodes(prediction.T,combined_tree.T)

    fitness_and_freqs = []
    for node in prediction.non_terminals:
        fitness_and_freqs.append([node.mean_fitness, node.polarizer]+list(node.mirror_node.temporal_frequency))

    fitness_and_freqs = np.array(fitness_and_freqs)
    return fitness_and_freqs

def clade_frequency_correlations_func_polarizer(params):
    # set up the prediction and pass all parameters to the wrapper function
    params.diffusion=1.0
    params.gamma = 1.0
    prediction = test_flu.predict_params(['polarizer'],
                                        params)


    # define the methodes for which the predictions are to be evaluated
    methods = [('polarizer', '_ext', prediction.terminals),
                ('polarizer', '_int', prediction.non_terminals)]
    distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)


    tbins_eval = [date(year=params.year, month = 6, day=1), date(year=params.year+1, month = 6, day=1)]
    combined_data = test_flu.make_combined_data(prediction, test_data, collapse = params.collapse)
    combined_tree = flu.flu_ranking(combined_data, time_bins = tbins_eval, pseudo_count = 0)
    combined_tree.expansion_score()

    tree_utils.find_internal_nodes(prediction.T,combined_tree.T)
    freqs = [node.mirror_node.temporal_frequency for node in prediction.non_terminals]
    polarizers= []
    for tau in mem_scale:
        prediction.calculate_polarizers(mem = tau)
        polarizers.append([node.polarizer for node in prediction.non_terminals])

    polarizers_and_freqs = np.hstack( (np.array(polarizers).T, np.array(freqs)))
    return polarizers_and_freqs


def clade_frequency_correlations_all_years(params, years = range(1995,2014)):
    fitness_and_freqs = []
    for year in years:
        params.year = year
        params.collapse = True
        params.sample_size =200
        params.omega=0.001
        fitness_and_freqs.append(clade_frequency_correlations_func_polarizer(params))
    return fitness_and_freqs

def plot_clade_frequencies_correlations(fitness_and_freqs, coli = 0):
    all_freqs = []
    plt.figure()
    for faf in fitness_and_freqs:
        tmp = faf[(faf[:,-2]<0.8)*(faf[:,-2]>0.005)]
        rank_fit = stats.rankdata(tmp[:,coli])
        rank_fit/=rank_fit.max()
        rank_freq = stats.rankdata(tmp[:,-1]/tmp[:,-2])
        rank_freq/=rank_freq.max()
        all_freqs.extend([(a,b) for a,b in zip(rank_fit,rank_freq)])
        plt.scatter(rank_fit, rank_freq)
        print stats.spearmanr(tmp[:,coli],tmp[:,-1]/tmp[:,-2])

    plt.xlabel('fitness rank')
    plt.ylabel('growth rank')
    all_freqs =  np.array(all_freqs)
    print 'all', stats.spearmanr(all_freqs[:,0], all_freqs[:,1])
    return all_freqs


if __name__=="__main__":

    # parse the commandline arguments
    parser = test_flu.make_flu_parser()
    params=parser.parse_args()
    if params.year>0:
        faf = clade_frequency_correlations_func_polarizer(params)
        plot_clade_frequencies_correlations([faf])
    else:
        coli=4
        faf = clade_frequency_correlations_all_years(params, years = range(2003,2014))
        normed_freqs = plot_clade_frequencies_correlations(faf, coli=coli)
        for ff in file_formats:
            plt.savefig(figure_folder+'fitness_frequency_correlations_years_2003_2014_gamma_'
                        +str(mem_scale[coli]/2)+ff)
        with open(analysis_folder+'fitness_frequency_correlations_years_2003_2014_gamma_'
                  +str(mem_scale[coli]/2)+'.pickle', 'w') as outfile:
            pickle.dump(normed_freqs, outfile)
        
        
