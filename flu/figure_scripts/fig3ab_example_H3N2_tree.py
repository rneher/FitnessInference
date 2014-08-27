#!/ebio/ag-neher/share/programs/bin/python2.7
#
#makes figure 3ab consisting of 2 trees of year 2007 colored with the prediction
#for the external nodes. It also produces a joint tree with the following season
#to illustrates the prediction. In the actual figure, the trees are rotated by 90degrees
#but without changing the scale. This script is based on analyze_HA1.py
#
import matplotlib
matplotlib.use('pdf')
import sys
sys.path.append('../../prediction_src')
import test_flu_prediction as flu
from Bio import Phylo,AlignIO,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os,argparse
from datetime import date

# set matplotlib plotting parameters
params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'lines.linewidth':1.5,
'text.usetex': True}
plt.rcParams.update(params)

# set flutype, prediction regions, and basic parameters
year=2007
flutype = 'H3N2'
prediction_regions = ["asia","north america"]
test_regions = ["asia","north america"]
sample_size = 50
boost = False
eps_branch_length = 1e-5
dscale = 5
D=0.5
collapse = False
aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'
#coding region of the HA1 sequence (the entire nucleotide sequence)
cds = {'begin':0, 'end':987, 'pad':0}

# make a name identifier that contains all the parameter information
base_name = flutype+'_'+str(year)+'_trees_pred_'+'_'.join(prediction_regions)+'_test_'+'_'.join(test_regions)
base_name=base_name.replace('north america', 'NA')
print base_name
name_mod = 'ssize_'+str(sample_size)+'_b_'+str(int(boost))+'_d_'+str(dscale)+'_D_'+str(D)
if collapse: name_mod+='_collapse'

# open annotations file
with open('../data/annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)
outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')
# load predictions by L&L to compare if desired
if os.path.isfile('../data/'+flutype+'_L_L_predictions.pickle'):
    with open('../data/'+flutype+'_L_L_predictions.pickle') as infile:
        laessig_prediction = pickle.load(infile)

#define time intervals use to predict and test
prediction_set={'start': date(year-1, 5,1), 'stop':date(year, 2,28),
                'regions':prediction_regions, 'sample_size':sample_size}
test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
            'regions':test_regions, 'sample_size':50}

# do the prediction
prediction = flu.predict(aln_fname, outgroup, annotation,\
                        ['mean_fitness'],
                        prediction_set, cds, time_bins = [], subsample_factor = 1.0, boost = boost,
                        eps_branch_length = 1e-5, collapse = collapse, dscale = dscale)

#define methods to evaluate the predictions and evaluate
methods = [('mean_fitness', '_ext', prediction.terminals),
            ('mean_fitness', '_int', prediction.non_terminals)]
distances, distances_epi, test_data = flu.evaluate(prediction, methods, test_set, aln_fname, outgroup, annotation,cds, flutype)
# add prediction if we have a useful one
if flutype=='H3N2' and (year in laessig_prediction) and (laessig_prediction[year].name not in [c.name for c in prediction.terminals]):
    otherseqsnames = [laessig_prediction[year].name]
else:
    otherseqsnames=[]
combined_data = flu.make_combined_data(prediction, test_data, otherseqsnames)

if flutype=='H3N2' and year in laessig_prediction:
    seq_labels = {prediction.best_node().name:'*'} #, laessig_prediction[year].name:"L&L"}
flu.label_nodes(prediction.T, seq_labels)
flu.label_nodes(combined_data.T, seq_labels)

# plot a combined figure
fig = plt.figure(figsize = (12,6))
#subplot 1: only the prediction data
ax = plt.subplot(121)
#add panel label
plt.text(-0.06,0.95,'A', transform = plt.gca().transAxes, fontsize = 36)
plt.title('until Feb '+str(year))
plt.tight_layout()
flu.plot_prediction_tree(prediction, axes=ax, cb=True, offset = 0.0005, internal=False)

#subplot 2: prediction data and test data
ax = plt.subplot(122)
#add panel label
plt.text(-0.06,0.95,'B', transform = plt.gca().transAxes, fontsize = 36)
plt.title('until Feb '+str(year)+' + season '+str(year)+"/"+str(year+1)+" (grey)")
pred_names = [c.name for c in prediction.terminals]
flu.erase_color(combined_data.T)
prediction.color_other_tree(combined_data.T.get_terminals(), offset = 0.0005)
for c in combined_data.T.get_terminals():
    if c.name in pred_names:
        pass
    #   c.color = (178, 34, 34 )
    else:
        c.color = (0,255,255)
        c.color = (100,100,100)
prediction.interpolate_color(combined_data.T)
flu.draw_flu_tree(combined_data.T, axes=ax, cb=False)

plt.savefig('../figures/Fig3AB_'+base_name+'_combined_'+name_mod+'.svg')


