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

# set up a toy data set using the classes in test_flu_prediction
prediction = test_flu.predict_params(['mean_fitness', 'ladder_rank', 'polarizer'], params)
base_name, name_mod = test_flu.get_fname(params)

methods = [ ('fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_int', prediction.non_terminals),
            ('ladder_rank', '', prediction.terminals),
            ('fitness','',prediction.terminals),
            ('polarizer', '_ext', prediction.terminals),
            ('polarizer', '_int', prediction.non_terminals)]  # last prediction is true fitness

distances, dist_epi_dummy, test_data = test_flu.evaluate(prediction, methods, params)
distances_list = [distances['average'],distances['minimal']]\
                    +[distances[m[0]+m[1]] for m in methods]

fitness_correlation_mf = stats.spearmanr([n.fitness for n in prediction.terminals],
                                      [n.mean_fitness for n in prediction.terminals])
fitness_correlation_pol = stats.spearmanr([n.fitness for n in prediction.terminals],
                                      [n.polarizer for n in prediction.terminals])
distances_list.extend(fitness_correlation_mf)
distances_list.extend(fitness_correlation_pol)
print distances_list

# write results to file
fname = '_'.join(map(str, [params.base_name, name_mod, 'prediction_results', 'dt', params.valdt*params.dt]))+'.dat'
methods = ['generation', 'average', 'minimal'] + [m[0]+m[1] for m in methods]
# if file does not exist, create and add header.
if not os.path.isfile(fname):
    with open(fname, 'w') as outfile:
        outfile.write('#'+'\t'.join(methods)+'\tcorrelation_mf\tpval_mf\tcorrelation_pol\tpval_pol\n')

with open(fname, 'a') as outfile:
    outfile.write('\t'.join(map(str,[a for a in [params.gen]+distances_list]))+'\n')

