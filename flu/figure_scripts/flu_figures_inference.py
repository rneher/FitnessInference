
#########################################################################################
#
# author: Richard Neher
# email: richard.neher@tuebingen.mpg.de
#
# Reference: Richard A. Neher, Colin A Russell, Boris I Shraiman. 
#            "Predicting evolution from the shape of genealogical trees"
#
##################################################
#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script that reads in precomputed repeated prediction of influenza and
#and plots the average predictions using external, internal nodes for each year 
#in addition, it compares this to predictions rewarding Koel et al mutations
#and to predictions using explicit temporal information (frequency dynamics within
#clades)
#
import glob,argparse,sys
sys.path.append('/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/src')
import test_flu_prediction as test_flu
import analysis_utils as AU
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
file_formats = [] #['.pdf', '.svg']

# set matplotlib plotting parameters
plt.rcParams.update(test_flu.mpl_params)
figure_folder = '../figures_ms/'

# set flutype, prediction regions, and basic parameters
parser = test_flu.make_flu_parser()
params=parser.parse_args()
params.year='????'
params.sample_size = 100
D = params.diffusion = 0.5
gamma = params.gamma = 3.0
omega = params.omega = 0.001
params.collapse = False
metric = 'nuc'

# make file identifiers
base_name, name_mod = test_flu.get_fname(params)
#remove year
base_name = '_'.join(base_name.split('_')[:1]+base_name.split('_')[2:])
base_name = base_name.replace('_????','')
# load data (with Koel boost and without), save in dictionary
prediction_distances={}
normed_distances={}
for boost in [0.0,0.5,1.0]:
    params.boost = boost
    years,tmp_pred, tmp_normed = AU.load_prediction_data(params, metric)
    prediction_distances.update(tmp_pred)
    normed_distances.update(tmp_normed)

##################################################################################
## main figure 3c
##################################################################################

# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal and L&L predictions
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    plt.plot([year-0.5, year+0.5], [prediction_distances[('minimal',boost,'minimal')][yi], 
                                    prediction_distances[('minimal',boost,'minimal')][yi]],
            lw=2, c='k', ls = '--')

for method, sym, col, shift, label  in [[('fitness,terminal nodes',0.0,'pred(T)'), 's', 'k', -0.25, 'top ranked terminal nodes'],
                                        [('fitness,internal nodes',0.0,'pred(I)'), 'd', 'r', 0.25, 'top ranked internal nodes ']]:
    plt.plot(years+shift, prediction_distances[method], sym, c= col, ms=8,
             label=label) #+r' $\bar{d}='+str(np.round(normed_distances[method][0],2))+'$')


# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel(r'$\Delta(\mathrm{prediction})$ to next season')
plt.xlabel('year')
plt.legend(loc=9, ncol=1,numpoints=1)
#add panel label
plt.text(-0.06,0.95,'C', transform = plt.gca().transAxes, fontsize = 36)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4C_'+base_name+'_'+name_mod+'_internal_external_revised'+ff)


##################################################################################
## Fig 4: compare bootstrap distributions of prediction results
## Bootstrapping is over years
##
##################################################################################

#sorted_methods = [a for a in sorted(normed_distances.items(), key=lambda x:x[1]) if a[0][0] 
#                  not in ['ladder rank', 'date', 'expansion, internal nodes', 'L&L'] or a[0][1]==0.0]
tick_labels = {    ('fitness,internal nodes', 0.0, 'pred(I)'):'internal',
                   ('fitness,terminal nodes', 0.0, 'pred(T)'):'terminal',
                   ('expansion, internal nodes', 0.0, 'growth'):'growth',
                   ('L&L', 0.0, r'L\&L'):r'L\&L',
                   ('ladder rank',0.0, 'ladder rank'):'ladder rank'}
sorted_methods = [a for a in sorted(normed_distances.items(), key=lambda x:x[1][0]) if a[0][:2] in 
                  [#('internal and expansion', 0.5),
                   #('internal and expansion', 0.0),
                   ('fitness,internal nodes', 0.0),
                   ('fitness,terminal nodes', 0.0),
                   ('expansion, internal nodes', 0.0),
                   ('L&L', 0.0),
                   ('ladder rank',0.0)] ]

plt.figure(figsize = (8,5))
plt.boxplot([a[1][1][-1] for a in sorted_methods],positions = range(len(sorted_methods)))
#plt.xticks(range(len(sorted_methods)), [a[0][-1] for a in sorted_methods], rotation=30, horizontalalignment='right')
plt.xticks(range(len(sorted_methods)), [tick_labels[a[0]] for a in sorted_methods], rotation=30, horizontalalignment='right')
plt.ylabel(r'distance $\bar{d}$ to next season')
plt.xlim([-0.5, len(sorted_methods)-0.5])
plt.grid()
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig5_'+base_name+'_'+name_mod+'_method_comparison'+ff)

##################################################################################
## Fig 3c-1 Comparison to L&L
##################################################################################
# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    plt.plot([year-0.5, year+0.5], [prediction_distances[('minimal',boost,'minimal')][yi],
                                    prediction_distances[('minimal',boost,'minimal')][yi]],
            lw=2, c='k', ls = '--')

method, sym, col, shift, label = ('fitness,terminal nodes',0.0,'pred(T)'), 's', 'k', -0.25, 'top ranked terminal nodes '
plt.plot(years+shift, prediction_distances[method], sym, c= col, ms=8, label=label+r' $\bar{d}='+str(np.round(normed_distances[method][0],2))+'$')
method, sym, col, shift, label = ('L&L',0.0,'L\&L'), 'o', 'r', 0.25, r'prediction by \L{}uksza and L\"assig'
plt.plot(years[AU.laessig_years(years)]+shift, prediction_distances[method][AU.laessig_years(years)],
         sym, c= col, ms=8, label=label+r' $\bar{d}='+str(np.round(normed_distances[method][0],2))+'$')

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
plt.text(0.02,0.9,'Fig.~3-S1', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4C_s1_'+base_name+'_'+name_mod+'_LL_external_revised'+ff)


##################################################################################
## Fig 3c-2 inclusion of Koel boost -- no temporal compnent
##################################################################################
# make figure
plt.figure(figsize = (12,6))
plt.title(r'Rewarding Koel mutations -- w/o calde growth estimate: $\bar{d}='
          +', '.join(map(str,[np.round(normed_distances[('fitness,internal nodes',boost,'pred(I)')][0],2)
                             for boost in [0.0, 0.5, 1.0]]))+'$ for $\delta = 0, 0.5, 1$', fontsize = 16)
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('fitness,internal nodes',0.0,'pred(I)'), 's', 'k', -0.25, 'pred(I)+Koel boost'
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    plt.plot([year-0.5, year+0.5], [prediction_distances[('minimal',boost,'minimal')][yi], 
                                    prediction_distances[('minimal',boost,'minimal')][yi]],
            lw=2, c='k', ls = '--')
    plt.plot(year+np.linspace(-0.5, 0.5,7)[1:-1:2], [prediction_distances[(method[0], koel, method[-1])][yi] for koel in [0.0, 0.5, 1.0]], 
         sym, c= col, ms=8,ls='-', label=label+r' $\bar{d}='+str(np.round(normed_distances[method][0],2))+'$')

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
    plt.savefig(figure_folder+'Fig4C_s2_'+base_name+'_'+name_mod+'_koel_boost_revised'+ff)

##################################################################################
## Fig 3c-3 inclusion of Koel boost -- with temporal compnent
##################################################################################
# make figure
plt.figure(figsize = (12,6))
plt.title(r'Rewarding Koel mutations -- with calde growth estimate: $\bar{d}='
          +', '.join(map(str,[np.round(normed_distances[('internal and expansion',boost,'pred(I)+growth')][0],2)
                             for boost in [0.0, 0.5, 1.0]]))+'$ for $\delta = 0, 0.5, 1$', fontsize = 16)

# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal 
method, sym, col, shift, label = ('internal and expansion',0.0,'pred(I)+growth'), 's', 'k', -0.25, 'pred(I)+Koel boost+ growth'
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
    plt.plot([year-0.5, year+0.5], [prediction_distances[('minimal',boost,'minimal')][yi], 
                                    prediction_distances[('minimal',boost,'minimal')][yi]],
            lw=2, c='k', ls = '--')
    plt.plot(year+np.linspace(-0.5, 0.5,7)[1:-1:2], [prediction_distances[(method[0], koel, method[-1])][yi] for koel in [0.0, 0.5, 1.0]], 
         sym, c= col, ms=8,ls='-', label=label+r' $\bar{d}='+str(np.round(normed_distances[method][0],2))+'$')

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
plt.text(0.02,0.93,'Fig.~3-S3', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4C_s3_'+base_name+'_'+name_mod+'_koel_boost_growth_revised'+ff)


##################################################################################
## Fig 3c-4 Temporal distribution of top strains
##################################################################################
# make figure
plt.figure(figsize = (10,6))
bins=np.cumsum([0,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31])
bc = (bins[1:]+bins[:-1])*0.5
bin_label = ['Jan', 'Feb', 'Mar', 'Apr', 'May','Jun','Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May']
params.boost = 0.0
sampling_dates = AU.load_date_distribution(params, 'mean_fitness')

for year in sorted(sampling_dates.keys()):
    y,x = np.histogram(sampling_dates[year], bins=bins)
    plt.plot(bc, y, 'o', label = str(year), ls='-')

# set limits, ticks, legends
plt.ylabel('distribution of predicted strains')
plt.xticks(bc, bin_label)
plt.xlabel('sampling date')
plt.legend(loc=1, ncol=2,numpoints=1)
#add panel label
plt.text(0.02,0.9,'Fig.~3-S4', transform = plt.gca().transAxes, fontsize = 20)
#save figure
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig3C_s4_'+base_name+'_'+name_mod+'_sampling_dates_by_year'+ff)

plt.figure()
all_dates = []
for d in sampling_dates.values(): all_dates.extend(d)
plt.hist(all_dates, bins=bins)
#add panel label
plt.text(0.02,0.9,'Fig.~3-S5', transform = plt.gca().transAxes, fontsize = 20)
plt.ylabel('distribution of predicted strains')
plt.xlabel('sampling date')
plt.xticks(bc, bin_label)
plt.tight_layout()
for ff in file_formats:
    plt.savefig(figure_folder+'Fig4C_s5_'+base_name+'_'+name_mod+'_sampling_dates'+ff)



