#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script provdies functions used to load and process the flu prediction data
#
import numpy as np
from scipy.stats import scoreatpercentile
import glob

# exclude some years from the comparison to the L&L predictions
# 2002 (illegitimate prediction of mislabled strain) and 2012 (no prediction availabe)
def laessig_years(years):
    return np.where((years!=2002)&(years!=2003)&(years<2012))[0]


methods = [(3,"fitness,terminal nodes", "pred(T)"), (4, "fitness,internal nodes", "pred(I)"), (5, "expansion, internal nodes", "growth"), 
           (6,"internal and expansion","pred(I)+growth"), (7,"external and time", "pred(T)+date"), (8,"ladder rank", "ladder rank"),(9,"date", "date") ]

def boot_strap(x,n=100):
    '''
    take an array and calculate the distribution of its mean by bootstrapping
    arguments:
        x   --  array whose elements are bootstrapped
        n   --  number of bootstrap replicates
    returns mean, (25th, median, 75th)precentile, boot_strap replicates
    '''
    bs = np.array([x[np.random.randint(len(x), size=len(x))].mean() for ni in xrange(n)])
    return bs.mean(), scoreatpercentile(bs,25), scoreatpercentile(bs, 50), scoreatpercentile(bs,75), bs

def load_prediction_data(flutype = 'H3N2', analysis_folder = '../analysis_may_feb', D=0.5, boost = 0.0,
              dscale = 5.0, prediction_regions = ["asia","north america"],
              test_regions = ["asia","north america"], sample_size = 100, metric = 'nuc'):

    # construct file mask and load files
    name_mod = 'ssize_'+str(sample_size)+'_b_'+str(boost)+'_d_'+str(dscale)+'_D_'+str(D)
    fmask = analysis_folder + flutype+'_year_????'+'_pred_'+'-'.join(prediction_regions)\
            +'_test_'+'-'.join(test_regions)+'_'+name_mod+'_'+metric+'.dat'
    fmask = fmask.replace(' ', '')
    flist = glob.glob(fmask)
    predictions = []
    for fname in flist:
        year = int(fname.split('_')[len(analysis_folder.split('_'))+1])
        predictions.append((year, np.loadtxt(fname)))

    predictions.sort(key = lambda x:x[0])
    # make a list of all years for which we have prediction
    years = np.array([a[0] for a in predictions])

    # calculate averages over the 50 replicates done for each year and normalize to random picks
    average_distance = np.array([a[1][:,0].mean()       for a in predictions])
    minimal_distance = np.array([a[1][:,1].mean()       for a in predictions])/average_distance
    prediction_distances = {(method,boost,label):np.array([a[1][:,methodi].mean()
                            for a in predictions])/average_distance
                            for methodi, method,label in methods}
    prediction_distances[('average',boost,'average')] = average_distance
    prediction_distances[('minimal',boost,'minimal')] = minimal_distance

    # normalize the observed distances to 
    normed_distances = {method:(np.mean((preds - minimal_distance)/(1-minimal_distance)),
                                boot_strap((preds - minimal_distance)/(1-minimal_distance),1000))
                       for method,preds in prediction_distances.iteritems()}


    # the L&L predictions sit in column 2
    prediction_distances[('L&L',boost,'L\&L')] = np.array([a[1][:,2].mean() for a in predictions])/average_distance

    normed_distances[('L&L',boost,'L\&L')] = (np.mean(((prediction_distances[('L&L', boost,'L\&L')] - minimal_distance)/(1-minimal_distance))[laessig_years(years)]),
                                      boot_strap(((prediction_distances[('L&L', boost,'L\&L')] - minimal_distance)/(1-minimal_distance))[laessig_years(years)])) 

    return years, prediction_distances, normed_distances

def load_date_distribution(flutype = 'H3N2', analysis_folder = '../analysis_may_feb', D=0.5, boost = 0.0,
              dscale = 5.0, prediction_regions = ["asia","north america"],
              test_regions = ["asia","north america"], sample_size = 100):
    '''
    returns the sampling dates of all predicted strains measured in days
    relative to Jan 1st of the year preceding the prediction year
    '''

    from datetime import date
    # construct file mask and load files
    name_mod = 'ssize_'+str(sample_size)+'_b_'+str(boost)+'_d_'+str(dscale)+'_D_'+str(D)
    fmask = analysis_folder + flutype+'_year_????'+'_topstrains_'+'-'.join(prediction_regions)\
            +'_test_'+'-'.join(test_regions)+'_'+name_mod+'.dat'
    fmask = fmask.replace(' ', '')
    flist = glob.glob(fmask)
    sampling_dates = {}
    for fname in flist:
        tmp_dates = []
        year = int(fname.split('_')[len(analysis_folder.split('_'))+1])
        base_line = date(year-1,1,1).toordinal()
        with open(fname, 'r') as infile:
            for line in infile:
                strain_year, strain_month, strain_day = map(int, line.split()[-1].split('-'))
                tmp_dates.append(date(strain_year, strain_month,strain_day).toordinal()-base_line)
        sampling_dates[year] = tmp_dates
    return sampling_dates
    
