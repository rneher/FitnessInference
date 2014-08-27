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
import test_flu_prediction as flu
from datetime import date
from Bio import Phylo, SeqIO
from matplotlib import pyplot as plt
import numpy as np
import glob, pickle
from scipy import stats

parser = argparse.ArgumentParser(description="test prediction algorithm on toy date")
parser.add_argument('--base_name', default='./', type=str, help='path to date including prefix')
parser.add_argument('--gen', default=2000, type=int, help='simulated generation used to center prediction interval')
parser.add_argument('--dt', default=100, type=int, help='number of simulated generations to sample sequences from')
parser.add_argument('--valdt', default=1, type=int, help='time interval to evalutation -- multiple of dt.')
parser.add_argument('--sample_size', default=100, type=int, help='max number of sequences to use')
parser.add_argument('--diffusion', default=0.5, type=float, help='Fitness diffusion coefficient')
parser.add_argument('--dscale', default=5.0, type=float, help='scale factor for time scale')
params = parser.parse_args()

# set up a toy data set using the classes in test_flu_prediction
cds = None
aln_fname = params.base_name+'_AT.fasta.gz'
prediction_regions = ['toy']
test_regions = ['toy']
prediction_set={}
test_set = {}
flutype = ''
ndt = params.valdt
with open(params.base_name+'_annotation.pickle', 'r') as infile:
    annotation = pickle.load(infile)

with open(params.base_name+'_outgroup.fa', 'r') as infile:
    outgroup = SeqIO.read(infile,'fasta')

prediction_set['start'] = params.gen-0.5*params.dt
prediction_set['stop'] =  params.gen+0.5*params.dt
prediction_set['regions'] = prediction_regions
prediction_set['sample_size'] = params.sample_size
test_set['start']= params.gen+(ndt - 0.25)*params.dt
test_set['stop'] = params.gen+(ndt + 0.25)*params.dt
test_set['regions'] = test_regions
test_set['sample_size'] = 50

prediction = flu.predict(aln_fname, outgroup, annotation,\
                        ['mean_fitness', 'ladder_rank'],
                        prediction_set, cds, time_bins = None, subsample_factor = 1.0, boost = 0.0,
                        eps_branch_length = 1e-5, collapse = False, dscale = params.dscale, D=params.diffusion)
methods = [('mean_fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_int', prediction.non_terminals),
            ('ladder_rank', '', prediction.terminals),
            ('fitness','',prediction.terminals)]  # last prediction is true fitness

distances, dist_epi_dummy, test_data = flu.evaluate(prediction, methods, test_set, aln_fname, outgroup,
                                                    annotation,cds, flutype=flutype)


mean_fitness_dist = test_data.mean_distance_to_sequence(prediction.best_node(method='mean_fitness').seq)
mean_fitness_int_dist = test_data.mean_distance_to_sequence(prediction.best_node(method='mean_fitness',
                                                            nodes = prediction.non_terminals).seq)
true_fitness_dist = test_data.mean_distance_to_sequence(prediction.best_node(method='fitness').seq)
branch_length_dist = test_data.mean_distance_to_sequence(prediction.best_node(method='ladder_rank').seq)
fitness_correlation = stats.spearmanr([n.fitness for n in prediction.terminals],
                                           [n.mean_fitness for n in prediction.terminals])

prediction_results = (params.gen, distances['average'], distances['minimal'], distances['fitness'],
                    distances['mean_fitness_ext'],distances['mean_fitness_int'],distances['ladder_rank'],
                    fitness_correlation[0], fitness_correlation[1])
print prediction_results[-1]

fname = params.base_name+'_prediction_results_'+str(params.sample_size)+'_d_'+str(params.dscale)\
    +'_D_'+str(params.diffusion)+'_dt_'+str(ndt*params.dt)+'.dat'
methods = ['generation', 'average', 'minimal','fitness_true', 'mean_fitness','mean_fitness_int', 'ladder_rank']
if not os.path.isfile(fname):
    with open(fname, 'w') as outfile:
        outfile.write('#'+'\t'.join(methods)+'\tcorrelation\tpval\n')

with open(fname, 'a') as outfile:
    outfile.write('\t'.join(map(str,[a for a in prediction_results]))+'\n')

