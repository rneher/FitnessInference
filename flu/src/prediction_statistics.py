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
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date

analysis_folder = test_flu.flu_analysis_folder
# parse the commandline arguments
parser = test_flu.make_flu_parser()
params=parser.parse_args()
params.pred = params.pred.replace('^',' ')
params.test =params.test.replace('^',' ')
params.subsample=0.7

# get run specific file names
fname_base, name_mod = test_flu.get_fname(params)

top_strain_method = 'mean_fitness'
# allocate arrays to save the predictions
nuc_dist_array = np.zeros((params.nreps, 12))
epi_dist_array = np.zeros((params.nreps, 12))
top_strains = []
for ii in xrange(params.nreps):
    # set up the prediction and pass all parameters to the wrapper function
    prediction = test_flu.predict_params(['mean_fitness', 'expansion_score', 'depth', 'polarizer',
                                          flu.combined_ranking_internal,
                                          flu.combined_ranking_external],
                                        params)

    # define the methodes for which the predictions are to be evaluated
    methods = [('mean_fitness', '_ext', prediction.terminals),
                ('mean_fitness', '_int', prediction.non_terminals),
                ('expansion_score', '_int', prediction.non_terminals),
                ('expansion_fitness', '', prediction.non_terminals),
                ('time_fitness', '', prediction.terminals),
                ('ladder_rank', '', prediction.terminals),
                ('date', '', prediction.terminals),
                ('polarizer', '_ext', prediction.terminals),
                ('polarizer', '_int', prediction.non_terminals)]
    distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)
    nuc_dist_array[ii,:] = [distances['average'],distances['minimal'],distances['L&L']]\
                    +[distances[m[0]+m[1]] for m in methods]
    epi_dist_array[ii,:] = [distances_epi['average'],distances_epi['minimal'],distances_epi['L&L']]\
                    +[distances_epi[m[0]+m[1]] for m in methods]
    # memorize the strain predicted best
    top_strains.append(prediction.best_node(method = top_strain_method, nodes = prediction.terminals))

#if file does not exist, create and write header
fname_nuc = analysis_folder+'_'.join([fname_base, name_mod, 'nuc.dat'])
if not os.path.isfile(fname_nuc):
    with open(fname_nuc, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_nuc, 'a') as outfile:
    np.savetxt(outfile, nuc_dist_array)

#if file does not exist, create and write header
fname_epi =analysis_folder+'_'.join([fname_base, name_mod, 'epi.dat'])
if not os.path.isfile(fname_epi):
    with open(fname_epi, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_epi, 'a') as outfile:
    np.savetxt(outfile, epi_dist_array)


# write the best nodes to file:
fname_strains = analysis_folder+'_'.join([fname_base, name_mod, top_strain_method, 'topstrains.dat'])
with open(fname_strains, 'a') as outfile:
    for strain in top_strains:
        outfile.write(strain.name+'\t'+str(strain.date)+'\n')

