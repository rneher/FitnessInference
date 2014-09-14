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

file_formats = ['.svg', '.pdf']
# set matplotlib plotting parameters
plt.rcParams.update(test_flu.mpl_params)

figure_folder = '../figures_ms/'
analysis_folder = test_flu.flu_analysis_folder
# parse the commandline arguments
parser = test_flu.make_flu_parser()
params=parser.parse_args()
params.year='????'
params.sample_size = 100
params.collapse = True

base_name, name_mod = test_flu.get_fname(params)
#remove year
base_name = '_'.join(base_name.split('_')[:1]+base_name.split('_')[2:])
base_name = base_name.replace('_????','')

# load data for different diffusion constants, distance_scales, and koel boosts
metric = 'nuc'
# divide by 2 to change units to pairwise distance
mem_time_scale = m_list = 2.0**np.arange(-6,4, 1) * 0.5
years, prediction_distance, normed_distance = AU.load_polarizer_data(params, metric)
LLii = AU.laessig_years(years)
tau_i=4
ti_ext = tau_i+3
ti_int = tau_i+3+len(mem_time_scale)
tau = mem_time_scale[tau_i]

##################################################################################
## Fig 5-1 polarizer figure
##################################################################################
# make figure showing the dependence on the scale parameter
fig = plt.figure(figsize = (12,9))
ax= plt.subplot(111)
plt.plot(mem_time_scale, np.ones_like(mem_time_scale)*normed_distance[LLii,1].mean(),
         c='k', lw=2, label = r"L\&L")

# plot values for terminal nodes, they reside in columns 2:(len(mem_time_scale)+2
plt.plot(mem_time_scale, np.ones_like(mem_time_scale)*normed_distance[:,2:(len(mem_time_scale)+2)].mean(axis=0),
         marker= 'o',ms=10, lw= 2, label = r'polarizer terminal')

# plot values for internal nodes, they reside in columns (len(mem_time_scale)+2):(2*len(mem_time_scale)+2
plt.plot(mem_time_scale, np.ones_like(mem_time_scale)*normed_distance[:,(len(mem_time_scale)+2):(2*len(mem_time_scale)+2)].mean(axis=0),
         marker= 'o',ms=10, lw= 2, label = r'polarizer internal')

plt.xscale('log')
ax.set_xlabel(r'time scale of memory [average pairwise distance]')
ax.set_ylabel(r'normalized distance $\bar{d}$')
plt.text(0.02,0.93,'Fig.~5-S1', transform = plt.gca().transAxes, fontsize = 20)
ax.set_ylim([0,1])
ax.set_xlim([min(mem_time_scale)*0.9,max(mem_time_scale)*1.1])
#plt.xscale('log')
plt.legend()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig5_s1_time_scale_dependence_polarizer_'+metric+ff)

##################################################################################
## Fig 4 polarizer figure -- individual years
##################################################################################
# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('polarizer,internal',0.0,''), 's', 'k', -0.25, 'pred(T)'
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    min_val = prediction_distance[yi,1]
    plt.plot([year-0.5, year+0.5], [min_val, min_val], lw=2, c='k', ls = '--')
    if yi: label_str = None, None
    else: label_str = 'terminal nodes '+r'$\tau='+str(tau)+'$', 'internal nodes '+r'$\tau='+str(tau)+'$'
    plt.plot(year-0.25,prediction_distance[yi,ti_ext], marker='o',
         c= 'r', ms=8,ls='', lw=2, label = label_str[0])
    plt.plot(year+0.25,prediction_distance[yi,ti_int] , marker = 'd',
         c= 'k', ms=8,ls='', lw=2, label = label_str[1])

# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel(r'$\Delta(\mathrm{prediction})$ to next season')
#plt.ylabel('nucleodide distance to next season\n(relative to average)')
plt.xlabel('year')
plt.legend(loc=9, ncol=1,numpoints=1)
#add panel label
plt.text(-0.06,0.95,'C', transform = plt.gca().transAxes, fontsize = 36)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4C_'+base_name+'_tau_'+str(tau)+'_polarizer_revised'+ff)


##################################################################################
## Fig 4-2 L&L comparison with polarizer
##################################################################################
# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('polarizer,internal',0.0,''), 's', 'k', -0.25, 'pred(T)'
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    min_val = prediction_distance[yi,1]
    plt.plot([year-0.5, year+0.5], [min_val, min_val], lw=2, c='k', ls = '--')
    if yi: label_str = None
    else: label_str = 'polarizer on terminal nodes '+r'$\tau='+str(tau)+'$'
    plt.plot(year-0.25,prediction_distance[yi,ti_ext], marker='o',
         c= 'r', ms=8,ls='', lw=2, label = label_str)

plt.plot(years[AU.laessig_years(years)]+0.25, prediction_distance[AU.laessig_years(years),2],
         marker = 'd', ls='',c= 'k', ms=8, label= r'\L{}uksza and L\"assig')

# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel(r'$\Delta(\mathrm{prediction})$ to next season')
#plt.ylabel('nucleodide distance to next season\n(relative to average)')
plt.xlabel('year')
plt.legend(loc=9, ncol=1,numpoints=1)
#add panel label
plt.text(0.02,0.93,'Fig.~4-S1', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4_s1_LL_'+base_name+'_tau_'+str(tau)+'_polarizer_revised'+ff)


##################################################################################
## Fig 4-3 varying time scale for polarizer
##################################################################################
# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('polarizer,internal',0.0,''), 's', 'k', -0.25, 'pred(T)'
terminal_indices = range(4, len(mem_time_scale)+3)
internal_indices = range(4+len(mem_time_scale),2*len(mem_time_scale)+3)
x_shifts = np.linspace(-0.4, 0.4,len(terminal_indices))
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    min_val = prediction_distance[yi,1]
    plt.plot([year-0.5, year+0.5], [min_val, min_val], lw=2, c='k', ls = '--')
    if yi: label_str = None, None
    else:
        trange = r'$\tau=2^{-6},\ldots,'+str(int(mem_time_scale[-1]))+'$'
        label_str = 'terminal nodes '+trange, 'internal nodes '+trange
    plt.plot(year+x_shifts,prediction_distance[yi,terminal_indices] ,
         c= 'r', ms=8,ls='-', lw=2, label = label_str[0])
    plt.plot(year+x_shifts,prediction_distance[yi,internal_indices] ,
         c= 'k', ms=8,ls='-', lw=2, label = label_str[1])

# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel(r'$\Delta(\mathrm{prediction})$ to next season')
#plt.ylabel('nucleodide distance to next season\n(relative to average)')
plt.xlabel('year')
plt.legend(loc=9, ncol=1,numpoints=1)
#add panel label
plt.text(0.02,0.93,'Fig.~4-S2', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4_s2_'+base_name+'_polarizer_revised'+ff)



##################################################################################
## Fig 5: compare bootstrap distributions of prediction results
## Bootstrapping is over years
##
##################################################################################

# load fitness prediction data
params.boost = 0.0
params.gamma = 3.0
params.omega = 0.1
params.diffusion = 0.5
params.collapse=False
years_I,prediction_distances_I, normed_distances_I = AU.load_prediction_data(params, metric)

plotted_methods =  {m:normed_distances_I[m] for m in [('expansion, internal nodes', 0.0, 'growth'),
                                                    ('L&L', 0.0, r'L\&L'),
                                                    ('ladder rank',0.0, 'ladder rank')]}

ti_ext_normed = tau_i+2
ti_int_normed = tau_i+2+len(mem_time_scale)
plotted_methods.update({('polarizer',tau,'external'):(normed_distance[:,ti_ext_normed].mean(), AU.boot_strap(normed_distance[:,ti_ext_normed], n=1000)),
                        ('polarizer',tau,'internal'):(normed_distance[:,ti_int_normed].mean(),AU.boot_strap(normed_distance[:,ti_int_normed], n=1000))})
tick_labels = {    ('fitness,internal nodes', 0.0, 'pred(I)'):'internal',
                   ('fitness,terminal nodes', 0.0, 'pred(T)'):'terminal',
                   ('expansion, internal nodes', 0.0, 'growth'):'growth',
                   ('L&L', 0.0, r'L\&L'):r'L\&L',
                   ('ladder rank',0.0, 'ladder rank'):'ladder rank',
                   ('polarizer',tau,'external'):r'terminal $\tau='+str(tau)+'$',
                   ('polarizer',tau,'internal'):r'internal $\tau='+str(tau)+'$'}
sorted_methods = [a for a in sorted(plotted_methods.items(), key=lambda x:x[1][0])]

plt.figure(figsize = (8,5))
plt.boxplot([a[1][1][-1] for a in sorted_methods],positions = range(len(sorted_methods)))
#plt.xticks(range(len(sorted_methods)), [a[0][-1] for a in sorted_methods], rotation=30, horizontalalignment='right')
plt.xticks(range(len(sorted_methods)), [tick_labels[a[0]] for a in sorted_methods], rotation=30, horizontalalignment='right')
plt.ylabel(r'distance $\bar{d}$ to next season')
plt.xlim([-0.5, len(sorted_methods)-0.5])
plt.grid()
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig5_'+base_name+'_polarizer_method_comparison'+ff)

   