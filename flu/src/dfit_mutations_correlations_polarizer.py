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
#script that analyzes a set of aligned sequences of the HA1 domain of influenza
#and attempts to predicts a strain close the population of the following year.
#We operate on historical data, hence the prediction is tested right away on
#next years strain collection. Optionally, it plots the reconstructed trees
#colored with the prediction ranking. In addition, it calculates delta fitness
#and the different types of mutations (syn, non-syn, Koel) for each branch and 
#writes them to a file.
#
################################################################################

import matplotlib
matplotlib.use('pdf')
import sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
import test_flu_prediction as test_flu
import predict_flu as flu
import tree_utils
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date

plt.rcParams.update(test_flu.mpl_params)

tree_figure_folder = '../figures_trees/'
analysis_folder = test_flu.flu_analysis_folder
# parse the commandline arguments
parser = test_flu.make_flu_parser()
parser.add_argument('--tau', default = 1.0, type = float, help= 'memory time scale of the tree polarizer')
params=parser.parse_args()
# get name snippets to link output files to run parameters
base_name, name_mod = test_flu.get_fname(params)
params.gamma=1.0
params.diffusion=1.0

# set up the prediction and pass all parameters to the wrapper function
prediction = test_flu.predict_params(['polarizer'], params)
prediction.calculate_polarizers(params.tau)

# define the methodes for which the predictions are to be evaluated
methods = [ ('polarizer', '_ext', prediction.terminals),
            ('polarizer', '_int', prediction.non_terminals)]
distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)

# calculate the fitness differentials for each internal branch and associate with 
# different types of mutations that happen on these branches
dfit = []
for node in prediction.non_terminals:
    for child in node.clades:
        delta_fitness = child.polarizer - node.polarizer
        muts = child.mutations
        aa_muts = child.aa_mutations
        if child.branch_length<0.01:
            dfit.append((delta_fitness,child.is_terminal(), len(muts), len(aa_muts),
                         len([m[0] for m in aa_muts if m[0] in flu.HA1_antigenic_sites]),
                         len([m[0] for m in aa_muts if m[0] in flu.cluster_positions])))
        else:
            print child,'excluded due to long branch length'
            print 'mutations:',muts
            print 'amino acids', aa_muts
dfit = np.array(dfit)
# write the array to file
fname = analysis_folder+base_name+'_tau_'+str(params.tau)+'_polarizer_dfit_mutations.dat'
with open(fname, 'w') as infile:
    infile.write('\t'.join(['#dfit', 'terminal', '#mutations', '#aa_muts', '#aa_muts at epitopes', '#aa_muts at Koel'])+'\n')
    np.savetxt(infile, dfit)
    

