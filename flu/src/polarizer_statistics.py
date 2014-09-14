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
params.diffusion = 1.0

# get run specific file names
fname_base, name_mod = test_flu.get_fname(params)

# allocate arrays to save the predictions
m_list = 2.0**np.arange(-6,4, 1)
nuc_dist_array = np.zeros((params.nreps, 3+2*len(m_list)))
epi_dist_array = np.zeros_like(nuc_dist_array)
top_strains = []
top_strain_method = 'polarizer'
for ii in xrange(params.nreps):
    # set up the prediction and pass all parameters to the wrapper function
    prediction = test_flu.predict_params(['polarizer'],params)
    test_data, test_set = test_flu.make_test_set(prediction, params)
    # define the methodes for which the predictions are to be evaluated
    methods = [ ('polarizer', '_ext', prediction.terminals),
                ('polarizer', '_int', prediction.non_terminals)]

    for mi,mem_time_scale in enumerate(m_list):
        prediction.calculate_polarizers(mem = mem_time_scale)
        distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params, test_data = test_data, test_set = test_set)
        nuc_dist_array[ii,mi+3] = distances['polarizer_ext']
        epi_dist_array[ii,mi+3] = distances_epi['polarizer_ext']
        nuc_dist_array[ii,mi+3+len(m_list)] = distances['polarizer_int']
        epi_dist_array[ii,mi+3+len(m_list)] = distances_epi['polarizer_int']
        
    nuc_dist_array[ii,:3] = [distances['average'],distances['minimal'],distances['L&L']]
    epi_dist_array[ii,:3] = [distances_epi['average'],distances_epi['minimal'],distances_epi['L&L']]
    # memorize the strain predicted best
    top_strains.append(prediction.best_node(method = top_strain_method, nodes = prediction.terminals))

#if file does not exist, create and write header
fname_nuc = analysis_folder+'_'.join([fname_base, 'polarizer_nuc.dat'])
if not os.path.isfile(fname_nuc):
    with open(fname_nuc, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_nuc, 'a') as outfile:
    np.savetxt(outfile, nuc_dist_array)

#if file does not exist, create and write header
fname_epi =analysis_folder+'_'.join([fname_base, 'polarizer_epi.dat'])
if not os.path.isfile(fname_epi):
    with open(fname_epi, 'w') as outfile:
        outfile.write('#average\tminimal\tL&L\t'+'\t'.join([m[0]+m[1] for m in methods])+'\n')
#append the results to existing file
with open(fname_epi, 'a') as outfile:
    np.savetxt(outfile, epi_dist_array)


# write the best nodes to file:
#fname_strains = analysis_folder+'_'.join([fname_base, name_mod, top_strain_method, 'topstrains.dat'])
#with open(fname_strains, 'a') as outfile:
#    for strain in top_strains:
#        outfile.write(strain.name+'\t'+str(strain.date)+'\n')

