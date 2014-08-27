#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script that reads in precomputed repeated prediction of influenza and
#and plots the average predictions using external, internal nodes for each year 
#in addition, it compares this to predictions rewarding Koel et al mutations
#and to predictions using explicit temporal information (frequency dynamics within
#clades)
#
import glob,argparse
import numpy as np
import matplotlib.pyplot as plt

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

# set flutype, prediction regions, and basic parameters
flutype = 'H3N2'
prediction_regions = ["asia","north america"]
test_regions = ["asia","north america"]
sample_size = 100
eps_branch_length = 1e-5
dscale = 5.0
D=0.5
collapse = False

# load data (with Koel boost and without), save in dictionary
predictions={}
for boost in [0,1]:
    name_mod = 'ssize_'+str(sample_size)+'_b_'+str(int(boost))+'_d_'+str(dscale)+'_D_'+str(D)
    if collapse: name_mod+='_collapse'
    fmask = '../analysis/' + flutype+'_year_????'+'_pred_'+'-'.join(prediction_regions)\
            +'_test_'+'-'.join(test_regions)+'_'+name_mod+'.dat'
    fmask = fmask.replace(' ', '')
    flist = glob.glob(fmask)
    predictions[boost] = []
    for fname in flist:
        year = int(fname.split('_')[2])
        predictions[boost].append((year, np.loadtxt(fname)))

    predictions[boost].sort(key = lambda x:x[0])

# make a list of all years for which we have prediction
years = np.array([a[0] for a in predictions[0]])
# exclude some years from the comparison to the L&L predictions
# 2002 (illegitimate prediction of mislabled strain) and 2012 (no prediction availabe)
laessig_years = np.where((years!=2002)&(years!=2012))[0]

# calculate averages over the 50 replicates done for each year and normalize to random picks
average_distance = np.array([a[1][:,0].mean()        for a in predictions[0]])
fit_ext_distance =    np.array([a[1][:,3].mean()  for a in predictions[0]])/average_distance
fit_int_distance =    np.array([a[1][:,4].mean()  for a in predictions[0]])/average_distance
boost_int_distance =    np.array([a[1][:,4].mean()  for a in predictions[1]])/average_distance
boost_comb_distance =    np.array([a[1][:,6].mean()  for a in predictions[1]])/average_distance
minimal_distance = np.array([a[1][:,1].mean()        for a in predictions[0]])/average_distance
laessig_distance =    np.array([a[1][:,2].mean()  for a in predictions[0]])/average_distance

# output prediction accuracy on a scale where optimal is 0 and random is 1
norm_distance_ext = np.round(np.mean((fit_ext_distance - minimal_distance)/(1-minimal_distance)),2)
norm_distance_int = np.round(np.mean((fit_int_distance- minimal_distance)/(1-minimal_distance)),2)
norm_distance_boost = np.round(np.mean((boost_int_distance - minimal_distance)/(1-minimal_distance)),2)
norm_distance_combi = np.round(np.mean((boost_comb_distance- minimal_distance)/(1-minimal_distance)),2)
norm_distance_LL = np.round(np.mean(((laessig_distance- minimal_distance)/(1-minimal_distance))[laessig_years]),2)
print "external:",      norm_distance_ext
print "internal:",      norm_distance_int
print "boosted:",       norm_distance_boost
print "combination:",   norm_distance_combi
print "L&L:",           norm_distance_LL

# make figure
plt.figure(figsize = (12,6))
# plot line for random expection
plt.plot([min(years)-0.5,max(years)+0.5], [1,1], lw=2, c='k')
# add shaded boxes and optimal and L&L predictions
for yi,year in enumerate(years):
    plt.gca().add_patch(plt.Rectangle([year-0.5, 0.2], 1.0, 1.8, color='k', alpha=0.05*(1+np.mod(year,2))))
#    if yi in laessig_years:
#        if yi==0:
#            plt.plot([year-0.5, year+0.5], [laessig_distance[yi], laessig_distance[yi]],
#                    lw=2, c='m',ls='-', label = r'\L{}uksza and L\"assig')
#        else:
#            plt.plot([year-0.5, year+0.5], [laessig_distance[yi], laessig_distance[yi]],
#                    lw=2, c='m',ls='-')
    plt.plot([year-0.5, year+0.5], [minimal_distance[yi], minimal_distance[yi]],
            lw=2, c='k', ls = '--')
plt.plot(years-0.0, fit_ext_distance,     's', c='k', ms=8, label='predicted distance')
#plt.plot(years-0.125, fit_int_distance,     '^', c='r', ms=8, label='predicted internal node $d='+str(norm_distance_int)+'$')
#plt.plot(years+0.125, boost_int_distance,   'd', c='b', ms=8, label='rewarding ``Koel" mutations $d='+str(norm_distance_boost)+'$')
#plt.plot(years+0.2, boost_comb_distance,  'v', c='r', ms=8, label='adding flu specifics $d='+str(norm_distance_combi)+'$')


# set limits, ticks, legends
plt.ylim([0.2, 1.7])
plt.yticks([0.5, 1, 1.5])
plt.xlim([min(years)-0.5,max(years)+0.5])
plt.xticks(years[::2])
plt.ylabel('distance relative to random pick')
plt.xlabel('year')
plt.legend(loc=9, ncol=1)
#add panel label
#plt.text(-0.06,0.95,'C', transform = plt.gca().transAxes, fontsize = 36)
#save figure
plt.tight_layout()
plt.savefig('../figures/Fig3C_simple_'+flutype+'_pred_'+'-'.join(prediction_regions).replace(' ', '')
            +'_comparison_ssize_'+str(sample_size)+'_d_'+str(dscale)+'_D_'+str(D)+'.svg')
