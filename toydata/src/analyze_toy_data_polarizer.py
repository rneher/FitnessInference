#!/ebio/ag-neher/share/programs/bin/python 
#!/usr/local/EPD/bin/python
#!/usr/bin/python
#########################################################################################
#
# author: Richard Neher
# email: richard.neher@tuebingen.mpg.de
#
# Reference: Richard A. Neher, Colin A Russell, Boris I Shraiman. 
#            "Predicting evolution from the shape of genealogical trees"
#
##################################################
# This scripts uses toy data as if it was flu data, i.e. uses similar data structures for
# annotation and the classes developed for flu to extract samples from the larger alignments.
##################################################

import sys,os, pickle, argparse
sys.path.append('../../flu/src')
import test_flu_prediction as test_flu
from datetime import date
from Bio import Phylo, SeqIO
from matplotlib import pyplot as plt
import numpy as np
import glob, pickle
from scipy import stats

parser = test_flu.make_toy_parser()
params = parser.parse_args()
params.pred = 'toy'
params.test = 'toy'
params.flutype='toy'
params.boost = 0.0
params.subsample = 1.0
params.year = params.gen
params.diffusion=1.0

# set up a toy data set using the classes in test_flu_prediction
m_list = 2.0**np.arange(-6,4, 1)
dist_array = np.zeros(2+2*len(m_list))
rho_array = np.zeros(len(m_list))

prediction = test_flu.predict_params(['polarizer'], params)
base_name, name_mod = test_flu.get_fname(params)
test_data, test_set = test_flu.make_test_set(prediction, params)

methods = [ ('polarizer', '_ext', prediction.terminals),
            ('polarizer', '_int', prediction.non_terminals)]

for mi,mem_time_scale in enumerate(m_list):
    prediction.calculate_polarizers(mem = mem_time_scale)
    distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params, test_data = test_data, test_set = test_set)
    dist_array[0]=distances['average']
    dist_array[1]=distances['minimal']
    dist_array[mi+2] = distances['polarizer_ext']
    dist_array[mi+2+len(m_list)] = distances['polarizer_int']
    rho_array[mi] = stats.spearmanr([n.fitness for n in prediction.terminals],
                                      [n.polarizer for n in prediction.terminals])[0]

# write results to file
fname = '_'.join(map(str, [params.base_name, 'prediction_results_polarizer', 'dt', params.valdt*params.dt]))+'.dat'
methods = ['generation', 'average', 'minimal'] + map(str, m_list) + map(str, m_list)
# if file does not exist, create and add header.
if not os.path.isfile(fname):
    with open(fname, 'w') as outfile:
        outfile.write('#'+'\t'.join(methods)+'\n')

with open(fname, 'a') as outfile:
    outfile.write('\t'.join(map(str,[a for a in [params.gen]+dist_array.tolist()]))+'\n')

# write results to file
fname = '_'.join(map(str, [params.base_name, 'corrcoeffs_polarizer', 'dt', params.valdt*params.dt]))+'.dat'
methods = ['generation'] + map(str, m_list)
# if file does not exist, create and add header.
if not os.path.isfile(fname):
    with open(fname, 'w') as outfile:
        outfile.write('#'+'\t'.join(methods)+'\n')

with open(fname, 'a') as outfile:
    outfile.write('\t'.join(map(str,[a for a in [params.gen]+rho_array.tolist()]))+'\n')

