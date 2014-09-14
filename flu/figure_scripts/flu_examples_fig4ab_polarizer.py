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
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/prediction_src')
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
# parse the commandline arguments
parser = test_flu.make_flu_parser()
parser.add_argument('--tau', default = 1.0, type = float, help= 'memory time scale of the tree polarizer')
params=parser.parse_args()
# get name snippets to link output files to run parameters
base_name, name_mod = test_flu.get_fname(params)
year = params.year
params.gamma=1.0
params.diffusion=1.0

# set up the prediction and pass all parameters to the wrapper function
prediction = test_flu.predict_params(['polarizer'],
                                    params)

prediction.calculate_polarizers(params.tau)
# define the methodes for which the predictions are to be evaluated
methods = [ ('polarizer', '_ext', prediction.terminals),
            ('polarizer', '_int', prediction.non_terminals)]
distances, distances_epi, test_data = test_flu.evaluate(prediction, methods, params)


# add prediction if we have a useful one
laessig_prediction = test_flu.get_LL(params.flutype)
if params.flutype=='H3N2' and (year in laessig_prediction) and (laessig_prediction[year].name not in [c.name for c in prediction.terminals]):
    otherseqsnames = [laessig_prediction[year].name]
else:
    otherseqsnames=[]
combined_data = test_flu.make_combined_data(prediction, test_data, otherseqsnames)

seq_labels = {prediction.best_node('polarizer').name:'*'} #, laessig_prediction[year].name:"L&L"}
tree_utils.label_nodes(prediction.T, seq_labels)
tree_utils.label_nodes(combined_data.T, seq_labels)

# plot a combined figure
fig = plt.figure(figsize = (12,6))
#subplot 1: only the prediction data
ax = plt.subplot(121)
#add panel label
plt.text(-0.06,0.95,'A', transform = plt.gca().transAxes, fontsize = 36)
plt.title('until Feb '+str(year))
plt.tight_layout()
tree_utils.plot_prediction_tree(prediction, method='polarizer', axes=ax, cb=True, offset = 0.0005, internal=False)

#subplot 2: prediction data and test data
ax = plt.subplot(122)
#add panel label
plt.text(-0.06,0.95,'B', transform = plt.gca().transAxes, fontsize = 36)
plt.title('until Feb '+str(year)+' + season '+str(year)+"/"+str(year+1)+" (grey)")
pred_names = [c.name for c in prediction.terminals]
tree_utils.erase_color(combined_data.T)
prediction.color_other_tree(combined_data.T.get_terminals(), method='polarizer', offset = 0.0005)
for c in combined_data.T.get_terminals():
    if c.name in pred_names:
        pass
    #   c.color = (178, 34, 34 )
    else:
        c.color = (0,255,255)
        c.color = (100,100,100)
prediction.interpolate_color(combined_data.T)
tree_utils.draw_tree(combined_data.T, axes=ax, cb=False)

for ff in file_formats:
    plt.savefig('../figures_ms/Fig4AB_'+base_name+'_combined_polarizer'+ff)


