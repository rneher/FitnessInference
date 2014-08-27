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

tree_figure_folder = '../figures_trees/'
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
parser.add_argument('--collapse', const = True, default=False, nargs='?', 
                    help='collapse internal branches with identical sequences')
parser.add_argument('--plot', const = True, default=False, nargs='?', help='plot trees')
parser.add_argument('--analysis', const = True, default=False, nargs='?', 
                    help='calculate delta fitness association tables and sampling histograms')


ssfactor = 1.0
params=parser.parse_args()
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
if params.collapse: name_mod+='_collapse'

# read in predictions by Luksza and Laessig
if os.path.isfile('../data/'+flutype+'_L_L_predictions.pickle'):
    with open('../data/'+flutype+'_L_L_predictions.pickle') as infile:
        laessig_prediction = pickle.load(infile)

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
    prediction_set={'start': date(year-1, 5,1), 'stop':date(year, 2,28),
                    'regions':prediction_regions, 'sample_size':params.sample_size}
    # test data start Oct 1st, end following year end of March
    test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                'regions':test_regions, 'sample_size':params.sample_size}
    
# define 6 month intervals to estimate changing clade frequencies.
tbins = [ date.fromordinal(prediction_set['stop'].toordinal()-ii*183) for ii in range(
        (prediction_set['stop'].toordinal()-prediction_set['start'].toordinal())//183,-1,-1)]

# set up the prediction and pass all parameters to the wrapper function
prediction = test_flu.predict(aln_fname, outgroup, annotation,\
                        ['mean_fitness', 'expansion_score', flu.combined_ranking_internal, flu.combined_ranking_external],
                        prediction_set, cds, time_bins = tbins, subsample_factor = ssfactor, boost = params.boost,
                        eps_branch_length = params.eps_branch, collapse = params.collapse, dscale = params.dscale, D=params.diffusion)

# define the methodes for which the predictions are to be evaluated
methods = [('mean_fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_int', prediction.non_terminals),
            ('expansion_score', '_int', prediction.non_terminals),
            ('expansion_fitness', '', prediction.non_terminals),
            ('time_fitness', '', prediction.terminals)]
distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, test_set, aln_fname, outgroup, annotation,cds, flutype)

# add the strain predicted by L&L if this happens to be a year they used
if flutype.startswith('H3N2') and (year in laessig_prediction) \
        and (laessig_prediction[year].name not in [c.name for c in prediction.terminals]):
    otherseqsnames = [laessig_prediction[year].name]
else:
    otherseqsnames=[]

# combine the test data, the prediction data and possible other sequences
# and build a tree
combined_data = test_flu.make_combined_data(prediction, test_data, otherseqsnames)

best_leaf = prediction.best_node()
print str(year)+'\t'+str(np.round((distances['mean_fitness_ext']-distances['minimal'])/(distances['average']-distances['minimal']),4)) + '\t' + best_leaf.name+'\n'

# make plots for leaves colored by inferred mean fitness
if params.plot:
    method = methods[0]
    internal = True
    if flutype.startswith('H3N2') and year in laessig_prediction:
        seq_labels = {prediction.best_node(method= method[0], nodes = method[2]).name:'*', laessig_prediction[year].name:"L&L"}
    else:
        seq_labels = {prediction.best_node().name:'*'}
        
    tree_utils.label_nodes(prediction.T, seq_labels)
    tree_utils.erase_color(prediction.T)
    tree_utils.erase_color(combined_data.T)
    tree_utils.plot_prediction_tree(prediction, internal=internal)
    plt.title("predicting "+flutype+" season "+str(year)+"/"+str(year+1)+": "+
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
    fname = analysis_folder+base_name+'_dfit_mutations_'+name_mod+'.dat'
    with open(fname, 'w') as infile:
        infile.write('\t'.join(['#dfit', 'terminal', '#mutations', '#aa_muts', '#aa_muts at epitopes', '#aa_muts at Koel'])+'\n')
        np.savetxt(infile, dfit)
        
    # determine the distributions of sampling dates and write to file    
    y,x,date_str = prediction.data.sampling_distribution()
    fname = analysis_folder+base_name+'_sampling_distribution_'+name_mod+'.pickle'
    with open(fname,'w') as outfile:
        pickle.dump((y,x,date_str), outfile)

