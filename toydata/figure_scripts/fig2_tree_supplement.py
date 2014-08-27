#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script that reads in precomputed analysis of simulated data and
#assemble a 3 panel figure illustrating and summarizing the performance
#of the prediction on simulated dat
#

import glob, sys,pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
sys.path.append('../../prediction_src')
import test_flu_prediction as flu
from Bio import Phylo

# set matplotlib plotting parameters
params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': False}
plt.rcParams.update(params)
# simulation data to use
data_path = '../data'
prefix = '20140418'
sigma=0.1

# indices of N,mutation, fliprate, sample_size, sampling interval, mean distance, minimal
# in the assempled files
ni = 0      # population size
mi = 1      # mutation rate
nflipi = 2  # flip rate at which position become advantageous
ssi = 3     # sample size
sdti = 4    # sampling interval
dscalei = 5 # stiffness parameter used in the analysis of the data
Di = 6 # stiffness parameter used in the analysis of the data
avgi = 9    # average distance to test data
mini = 10    # minimal distance to prediction
corri = 15  # rank correlation coefficient between true and estimated fitness

# make figures
N=10000
L=2000
nflips = .03
ssize =200
dscale = 3.0
D=0.5
mu=0.0001
sample_size = 200
# run example prediction
print "make tree plots"
plt.figure(figsize = (12,5))
# make example figure
######################################################
## instantaneous sampling
######################################################
#label positions
xpos = -0.3
ypos = 0.95
ax = plt.subplot(121)

sdt = 100
file_base = data_path+'_'.join(map(str,['/N', N, 'L', L, 'nflip',nflips,'mu',mu,'sdt', sdt]))+\
    '/'+prefix+'_seqs'
with open(file_base+'_annotation.pickle', 'r') as infile:
    annotation = pickle.load(infile)

with open(file_base+'_outgroup.fa', 'r') as infile:
    from Bio import SeqIO
    outgroup = SeqIO.read(infile,'fasta')

print file_base, outgroup
aln_fname = file_base+'_AT.fasta.gz'
cds = {'begin':0, 'end':0, 'pad':0}
dt=100
year = 6000 # 
#year = 6800
#year = 5500
#year = 9500

prediction_set={'start':year-0.5*dt, 'stop':year+0.5*dt, 'regions': ['toy'], 'sample_size':sample_size}
prediction1 = flu.predict(aln_fname, outgroup, annotation,\
                        ['mean_fitness'],
                        prediction_set, cds, time_bins = None, subsample_factor = 1.0, boost = 0.0,
                        eps_branch_length = 1e-5, collapse = False, dscale = dscale, D=D)
plt.text(xpos,ypos,'A', transform = plt.gca().transAxes, fontsize = 36)
Phylo.draw(prediction1.T, axes= ax, label_func= lambda x:'', show_confidence = False)


###########################################################
## continuous sampling
##########################################################
sdt = 1
file_base = data_path+'_'.join(map(str,['/N', N, 'L', L, 'nflip',nflips,'mu',mu,'sdt', sdt]))+\
    '/'+prefix+'_seqs'
with open(file_base+'_annotation.pickle', 'r') as infile:
    annotation = pickle.load(infile)

with open(file_base+'_outgroup.fa', 'r') as infile:
    from Bio import SeqIO
    outgroup = SeqIO.read(infile,'fasta')

print file_base, outgroup
aln_fname = file_base+'_AT.fasta.gz'
cds = {'begin':0, 'end':0, 'pad':0}

prediction2 = flu.predict(aln_fname, outgroup, annotation,\
                        ['mean_fitness'],
                        prediction_set, cds, time_bins = None, subsample_factor = 1.0, boost = 0.0,
                        eps_branch_length = 1e-5, collapse = False, dscale = dscale, D=D)

# make summary figure
ax = plt.subplot(122)
plt.text(xpos,ypos,'B', transform = plt.gca().transAxes, fontsize = 36)
Phylo.draw(prediction2.T, axes= ax, label_func= lambda x:'', show_confidence = False)

plt.tight_layout()
plt.savefig('../figures/Fig2_supp_tree_comparison.pdf')

