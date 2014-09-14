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
params=parser.parse_args()
# get name snippets to link output files to run parameters
base_name, name_mod = test_flu.get_fname(params)
# set up the prediction and pass all parameters to the wrapper function
prediction = test_flu.predict_params(['mean_fitness','expansion_score', 'depth', 'polarizer'],
                                    params)

# define the methodes for which the predictions are to be evaluated
methods = [('mean_fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_int', prediction.non_terminals),
            ('expansion_score', '_int', prediction.non_terminals),
            ('expansion_fitness', '', prediction.non_terminals),
            ('time_fitness', '', prediction.terminals),
            ('polarizer', '_ext', prediction.terminals),
            ('polarizer', '_int', prediction.non_terminals)]
distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)


# make plots for leaves colored by inferred mean fitness
if params.plot:
    laessig_prediction = test_flu.get_LL(params.flutype)
    year = params.year
    # add the strain predicted by L&L if this happens to be a year they used
    if params.flutype.startswith('H3N2') and (year in laessig_prediction) \
            and (laessig_prediction[year].name not in [c.name for c in prediction.terminals]):
        otherseqsnames = [laessig_prediction[year].name]
    else:
        otherseqsnames=[]
    
    # combine the test data, the prediction data and possible other sequences
    # and build a tree
    combined_data = test_flu.make_combined_data(prediction, test_data, otherseqsnames)

    method = methods[0]
    color_internal = True
    
    if params.flutype.startswith('H3N2') and year in laessig_prediction:
        seq_labels = {prediction.best_node(method= method[0], nodes = method[2]).name:'*', laessig_prediction[year].name:"L&L"}
    else:
        seq_labels = {prediction.best_node().name:'*'}
        
    tree_utils.label_nodes(prediction.T, seq_labels)
    tree_utils.erase_color(prediction.T)
    tree_utils.erase_color(combined_data.T)
    tree_utils.plot_prediction_tree(prediction, internal=color_internal)
    plt.title("predicting "+params.flutype+" season "+str(year)+"/"+str(year+1)+": "+
              str(np.round(distances[method[0]+method[1]]/distances['average'],4)))
    #plt.savefig('../figures/'+base_name+'_prediction_'+method[0]+method[1]+'_'+name_mod+'.pdf')
    
    # plot a combined figure
    # color according to season
    pred_names = [c.name for c in prediction.terminals]
    tree_utils.erase_color(combined_data.T)
    for c in combined_data.T.get_terminals():
        if c.name in pred_names:
            c.color = (178, 34, 34 )
        else:
            c.color = (0,255,255)
    prediction.interpolate_color(combined_data.T)
    
    tree_utils.label_nodes(combined_data.T, seq_labels)
    fig = plt.figure(figsize = (10,6))
    ax = plt.subplot(121)
    plt.title('until Feb '+str(year))
    tree_utils.draw_tree(prediction.T, axes=ax, cb=True)
    ax = plt.subplot(122)
    plt.title('before Mar '+str(year)+' (red) + season '+str(year)+"/"+str(year+1)+" (cyan)")
    tree_utils.draw_tree(combined_data.T, axes=ax, cb=False)
    plt.tight_layout()
    plt.savefig(tree_figure_folder+base_name+'_'+method[0]+method[1]+'_'+name_mod+'.pdf')


# perform additional analysis on the tree and the prediction
if params.analysis:
    # calculate the fitness differentials for each internal branch and associate with 
    # different types of mutations that happen on these branches
    dfit = []
    for node in prediction.non_terminals:
        for child in node.clades:
            delta_fitness = child.mean_fitness - node.mean_fitness
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
    fname = analysis_folder+base_name+'_'+name_mod+'_dfit_mutations.dat'
    with open(fname, 'w') as infile:
        infile.write('\t'.join(['#dfit', 'terminal', '#mutations', '#aa_muts', '#aa_muts at epitopes', '#aa_muts at Koel'])+'\n')
        np.savetxt(infile, dfit)
        
    # determine the distributions of sampling dates and write to file    
    y,x,date_str = prediction.data.sampling_distribution()
    fname = analysis_folder+base_name+'_'+name_mod+'_sampling_distribution.pickle'
    with open(fname,'w') as outfile:
        pickle.dump((y,x,date_str), outfile)

