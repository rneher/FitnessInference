#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script that reads in precomputed repeated prediction of influenza and
#and plots the average prediction quality as a function of the diffusion constant and the
#scale parameter gamma.
#
#
import glob,argparse
import numpy as np
import matplotlib.pyplot as plt
import analysis_utils as AU

file_formats = ['.svg', '.pdf']
# set matplotlib plotting parameters
params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
plt.rcParams.update(params)

figure_folder = '../figures_ms/'
analysis_folder = '../analysis_may_feb/'
# set flutype, prediction regions, and basic parameters
flutype = 'H3N2'
prediction_regions = ["asia","north america"]
test_regions = ["asia","north america"]
sample_size = 100
Dlist = [0.2, 0.5]
glist = [0.5,1.0,2.0,3.0, 5.0]

# load data for different diffusion constants, distance_scales, and koel boosts
prediction_distance={}
normed_distance={}
metric = 'nuc'
for D in Dlist:
    for dscale in glist:
        prediction_distance[(D,dscale)]={}
        normed_distance[(D,dscale)]={}
        for boost in [0.0,0.5,1.0]:
            years, tmp_pred, tmp_normed = AU.load_prediction_data(analysis_folder=analysis_folder, D=D,
                                                           dscale=dscale,boost=boost, sample_size=sample_size,
                                                           metric=metric)
            prediction_distance[(D,dscale)].update(tmp_pred)
            normed_distance[(D,dscale)].update(tmp_normed)

# make figure showing the dependence on the scale parameter
fig = plt.figure(figsize = (12,9))
ax= plt.subplot(111)
plt.plot(glist, np.ones_like(glist)*normed_distance[(0.5,5.0)][('L&L',0,'L\&L')][0], c='k', lw=2, label = r"L\&L")
boost = 0.0
for di,D in enumerate(Dlist):
    plt.plot(glist, [normed_distance[(D,dscale)][('fitness,terminal nodes',boost,'pred(T)')][0] for dscale in glist],
             marker= 'o',ms=10, lw= 2, label = 'top ranked terminal nodes $D='+str(D)+'$')

boost = 0.5
for di,D in enumerate(Dlist):
    plt.plot(glist, [normed_distance[(D,dscale)][('internal and expansion',boost,'pred(I)+growth')][0] for dscale in glist],
             marker= 'o',ms=10, lw= 2, label = 'Internal nodes + Koel('+str(boost)+') + growth $D='+str(D)+'$')

ax.set_xlabel(r'scale parameter $\gamma$')
ax.set_ylabel(r'normalized distance $\bar{d}$')
plt.text(0.02,0.93,'Fig.~4-S1', transform = plt.gca().transAxes, fontsize = 20)
ax.set_ylim([0,1])
ax.set_xlim([min(glist)-0.5,max(glist)+0.5])
plt.legend()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4_S1_D_gamma_dependence_'+metric+ff)
