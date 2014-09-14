#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script that reads in precomputed repeated prediction of influenza and
#and plots the average prediction quality as a function of the diffusion constant and the
#scale parameter gamma.
#
#
import glob,argparse,sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
import test_flu_prediction as test_flu
import numpy as np
import matplotlib.pyplot as plt
import analysis_utils as AU

file_formats = [] #['.svg', '.pdf']
# set matplotlib plotting parameters
plt.rcParams.update(test_flu.mpl_params)

figure_folder = '../figures_ms/'
analysis_folder = test_flu.flu_analysis_folder
# parse the commandline arguments
parser = test_flu.make_flu_parser()
params=parser.parse_args()
params.year='????'
params.sample_size = 100
Dlist = [0.2, 0.5]
glist = [1.0,2.0,3.0, 5.0]
olist = [0.1]
boost = 0.0

base_name, name_mod = test_flu.get_fname(params)
#remove year
base_name = '_'.join(base_name.split('_')[:1]+base_name.split('_')[2:])
base_name = base_name.replace('_????','')
LLy = AU.laessig_years(years)

params.collapse = False

# load data for different diffusion constants, distance_scales, and koel boosts
prediction_distance={}
normed_distance={}
metric = 'nuc'
for D in Dlist:
    for gamma in glist:
        for omega in olist:
            params.diffusion, params.gamma, params.omega = D,gamma, omega
            prediction_distance[(D,gamma,omega)]={}
            normed_distance[(D,gamma,omega)]={}
            years, tmp_pred, tmp_normed = AU.load_prediction_data(params, metric)
            prediction_distance[(D,gamma,omega)].update(tmp_pred)
            normed_distance[(D,gamma,omega)].update(tmp_normed)

# make figure showing the dependence on the scale parameter
fig = plt.figure(figsize = (12,9))
ax= plt.subplot(111)
plt.plot(glist, np.ones_like(glist)*normed_distance[(Dlist[0],glist[0],olist[0])][('L&L',boost,'L\&L')][0], c='k', lw=2, label = r"L\&L")

for omega in olist:
    for di,D in enumerate(Dlist):
        label_mod = r' $D='+str(D)+r'$'
        plt.plot(glist, [normed_distance[(D,gamma,omega)][('fitness,terminal nodes',boost,'pred(T)')][0] for gamma in glist],
                 marker= 'o',ms=10, lw= 2, label = r'top ranked terminal nodes'+label_mod)

 #       plt.plot(glist, [normed_distance[(D,gamma,omega)][('polarizer,terminal',boost,'')][0] for gamma in glist],
 #                marker= 'o',ms=10, lw= 2, ls=':', label = r'polarizer external'+label_mod)

        plt.plot(glist, [normed_distance[(D,gamma,omega)][('fitness,internal nodes',boost,'pred(I)')][0] for gamma in glist],
                 marker= 'o',ms=10, lw= 2, ls='--', label = r'top ranked internal nodes'+label_mod)

        plt.plot(glist, [normed_distance[(D,gamma,omega)][('expansion, internal nodes', 0.0, 'growth')][0] for gamma in glist],
                 marker= 'o',ms=10, lw= 2, ls='-.', label = r'expansion'+label_mod)
#        plt.plot(glist, [normed_distance[(D,gamma,omega)][('polarizer,internal',boost,'')][0] for gamma in glist],
#                 marker= 'o',ms=10, lw= 2, ls=':', label = r'polarizer interal'+label_mod)
#    boost = 0.0
#    for di,D in enumerate(Dlist):
#        plt.plot(glist, [normed_distance[(D,gamma,omega)][('internal and expansion',boost,'pred(I)+growth')][0] for gamma in glist],
#                 marker= 'o',ms=10, lw= 2, label = 'Internal nodes + Koel('+str(boost)+') + growth'+label_mod)

ax.set_xlabel(r'scale parameter $\gamma$')
ax.set_ylabel(r'normalized distance $\bar{d}$')
plt.text(0.02,0.93,'Fig.~4-S1', transform = plt.gca().transAxes, fontsize = 20)
ax.set_ylim([0,1])
ax.set_xlim([min(glist)*0.9,max(glist)*1.1])
#plt.xscale('log')
plt.legend()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4_S1_D_gamma_dependence_'+metric+ff)


##################################################################################
## Fig 4-2 varying gamma
##################################################################################
# make figure
plt.figure(figsize = (12,6))
boost=0.0
D=0.2
title_str = r'Varying $\gamma:\; \bar{d}='\
          +', '.join(map(str,[np.round(normed_distance[(D,gamma,omega)][('fitness,terminal nodes',boost,'pred(T)')][0],2)\
                             for gamma in glist]))+'$'  #+r' $ for $\gamma = 1.0, 2.0, 3.0, 5.0$'
#plt.title(title_str, fontsize = 16)
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('fitness,terminal nodes',0.0,'pred(T)'), 's', 'k', -0.25, 'pred(T)'
method, sym, col, shift, label = ('polarizer,internal',0.0,''), 's', 'k', -0.25, 'pred(T)'
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    plt.plot([year-0.5, year+0.5], [prediction_distance[(D,gamma,omega)][('minimal',boost,'minimal')][yi], 
                                    prediction_distance[(D,gamma,omega)][('minimal',boost,'minimal')][yi]],
            lw=2, c='k', ls = '--')
    plt.plot(year+np.linspace(-0.5, 0.5,9)[1:-1:2], [prediction_distance[(D,gamma,omega)][(method[0], boost, method[-1])][yi] for gamma in glist],
         sym, c= col, ms=8,ls='-', label=label+r' $\bar{d}='+str(np.round(normed_distance[(D,gamma,omega)][method][0],2))+'$')

# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel(r'$\Delta(\mathrm{prediction})$ to next season')
#plt.ylabel('nucleodide distance to next season\n(relative to average)')
plt.xlabel('year')
#plt.legend(loc=9, ncol=1,numpoints=1)
#add panel label
plt.text(0.02,0.93,'Fig.~3-S2', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4_s2_'+base_name+'_'+name_mod+'_gamma_revised'+ff)


