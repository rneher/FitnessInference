#########################################################################################
#
# author: Richard Neher
# email: richard.neher@tuebingen.mpg.de
#
# Reference: Richard A. Neher, Colin A Russell, Boris I Shraiman. 
#            "Predicting evolution from the shape of genealogical trees"
#
##################################################
# WRAPPER FUNCTIONS TO PREDICT AND EVALUATE THE PREDICTIONS
##################################################


import pickle,os, argparse
from datetime import date
import predict_flu as PF
import numpy as np
from Bio import SeqIO

# matplotlib params
mpl_params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}

flu_analysis_folder = '/ebio/ag-neher/share/users/rneher/FluPrediction_code/flu/analysis_may_feb/'
if not os.path.isdir(flu_analysis_folder):
    os.mkdir(flu_analysis_folder)

def make_flu_parser():
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
    parser.add_argument('--gamma', default=5.0, type=float, help='scale factor for time scale')
    parser.add_argument('--omega', default=0.001, type=float, help='sampling fraction')
    parser.add_argument('--collapse', const = True, default=False, nargs='?', 
                        help='collapse internal branches with identical sequences')
    parser.add_argument('--plot', const = True, default=False, nargs='?', help='plot trees')
    parser.add_argument('--analysis', const = True, default=False, nargs='?', 
                        help='calculate delta fitness association tables and sampling histograms')
    parser.add_argument('--subsample', default = 1.0, type = float, help='fraction of data to use if it exceeds sample size')
    parser.add_argument('--nreps', default='1', type=int, help='number of repeated runs using different subsamples of the data')
    return parser

def make_toy_parser():
    #parse the command line arguments for a toy data run
    parser = argparse.ArgumentParser(description="test prediction algorithm on toy date")
    parser.add_argument('--base_name', default='./', type=str, help='path to date including prefix')
    parser.add_argument('--gen', default=2000, type=int, help='simulated generation used to center prediction interval')
    parser.add_argument('--dt', default=100, type=int, help='number of simulated generations to sample sequences from')
    parser.add_argument('--valdt', default=1, type=int, help='time interval to evalutation -- multiple of dt.')
    parser.add_argument('--sample_size', default=100, type=int, help='max number of sequences to use')
    parser.add_argument('--diffusion', default=0.5, type=float, help='Fitness diffusion coefficient')
    parser.add_argument('--gamma', default=5.0, type=float, help='scale factor for time scale')
    parser.add_argument('--omega', default=0.001, type=float, help='sampling fraction')
    parser.add_argument('--eps_branch', default=1e-5, type=float, help='minimal branch length for inference')
    parser.add_argument('--collapse', const = True, default=False, nargs='?', 
                        help='collapse internal branches with identical sequences')
    return parser

def get_fname(params):
    '''
    parse the parameter name space and produce a data specific and an analysis specific modifier for
    to link output files to this parameter set
    '''
    prediction_regions = params.pred.split(',')
    test_regions = params.test.split(',')
    base_name = '_'.join(map(str,[params.flutype,'year',params.year, 'pred']+ 
                         ["-".join(prediction_regions)]+['test']+
                         ['-'.join(test_regions), 'ssize', params.sample_size])).replace(' ', '')
    name_mod = '_'.join(map(str, ['b', params.boost, 'd', params.gamma,
                                  'D',params.diffusion, 'w', params.omega]))
    if params.collapse: name_mod+='_collapse'    
    return base_name, name_mod


def get_LL(flutype):
    # read in predictions by Luksza and Laessig
    if os.path.isfile('../data/'+flutype+'_L_L_predictions.pickle'):
        with open('../data/'+flutype+'_L_L_predictions.pickle') as infile:
            laessig_prediction = pickle.load(infile)
    else:
        laessig_prediction=None
    return laessig_prediction

def get_aln_annotation_outgroup_cds(params):
    aln_fname = '../data/'+params.flutype+'_HA1_all_years_filtered.fasta.gz'
    
    if params.flutype.startswith('H3N2'):
        cds = {'begin':0, 'end':987, 'pad':0}
    else:
        cds = {'begin':0, 'end':300, 'pad':0}
        
    # open annotations file
    with open('../data/'+params.flutype+'_annotations.pickle', 'r') as infile:
        annotation = pickle.load(infile)
    outgroup = SeqIO.read('../data/'+params.flutype+'_outgroup.fasta', 'fasta')

    return aln_fname, annotation, outgroup, cds

def get_aln_annotation_outgroup_cds_toy(params):
    aln_fname = params.base_name+'_AT.fasta.gz'
    cds=None
    with open(params.base_name+'_annotation.pickle', 'r') as infile:
        annotation = pickle.load(infile)
    
    with open(params.base_name+'_outgroup.fa', 'r') as infile:
        outgroup = SeqIO.read(infile,'fasta')
        
    return aln_fname, annotation, outgroup, cds


def predict_params(methods, params):
    flutype = params.flutype
    if 'toy'!=params.pred:
        year=params.year
    prediction_regions = params.pred.split(',')
    test_regions = params.test.split(',')
    pseudo_count = 5
        
    # define prediction and test data sets by choosing start and stop dates 
    if "oceania" in prediction_regions:
        prediction_set={'start': date(year-2, 10,1), 'stop':date(year-1, 9,30),
                        'regions':prediction_regions, 'sample_size':params.sample_size}
    elif "toy" in prediction_regions:
        prediction_set={}
        prediction_set['start'] = params.gen-0.5*params.dt
        prediction_set['stop'] =  params.gen+0.5*params.dt
        prediction_set['regions'] = prediction_regions
        prediction_set['sample_size'] = params.sample_size
    else:  # assume northern hemisphere
        # begin previous year on may 1st, end this year on Feb 28
        prediction_set={'start': date(year-1,5,1), 'stop':date(year, 2,28),
                        'regions':prediction_regions, 'sample_size':params.sample_size}
    
    if 'toy' in params.pred:
        aln_fname, annotation, outgroup,cds = get_aln_annotation_outgroup_cds_toy(params)
        tbins = None
    else:
        aln_fname, annotation, outgroup,cds = get_aln_annotation_outgroup_cds(params)
        # define 105 day intervals to estimate changing clade frequencies.
        # chosen to have 3 intervals between May and Feb
        bin_dt = 105
        tbins = [ date.fromordinal(prediction_set['stop'].toordinal()-ii*bin_dt) for ii in range(
                (prediction_set['stop'].toordinal()-prediction_set['start'].toordinal())//bin_dt,-1,-1)]

    prediction_data = PF.flu_alignment(aln_fname,outgroup,annotation,cds=cds, subsample_factor = params.subsample,
                                      criteria = [[prediction_set['start'], prediction_set['stop'], 
                                    [reg],prediction_set['sample_size']] for reg in prediction_set['regions']], collapse=params.collapse)

    # PREDICT
    prediction = PF.flu_ranking(prediction_data, eps_branch_length = params.eps_branch, boost = params.boost,
                                time_bins = tbins, methods=methods, D=params.diffusion,samp_frac = params.omega,
                                 distance_scale = params.gamma,pseudo_count = pseudo_count)
    prediction.predict()
    print "prediction done"
    return prediction

def make_test_set(prediction, params):
    flutype = params.flutype
    if 'toy'!=params.test:
        year=params.year
    prediction_regions = params.pred.split(',')
    test_regions = params.test.split(',')
    # define test data sets by choosing start and stop dates 
    if "oceania" in test_regions:
        test_set = {'start':date(year, 3,1), 'stop':date(year, 9,30),
                    'regions':test_regions, 'sample_size':params.sample_size}
    elif 'toy' in params.test:
        test_set = {}
        test_set['start']= params.gen+(params.valdt - 0.25)*params.dt
        test_set['stop'] = params.gen+(params.valdt + 0.25)*params.dt
        test_set['regions'] = test_regions
        test_set['sample_size'] = 50        
    else:  # assume northern hemisphere
        test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                    'regions':test_regions, 'sample_size':params.sample_size}
    
    # load data to test
    test_data = PF.flu_alignment(prediction.data.aln_file_name,prediction.data.outgroup,
                                 prediction.data.annotation,cds=prediction.data.cds,
                             criteria = [[test_set['start'], test_set['stop'],
                            [reg],test_set['sample_size']] for reg in test_set['regions']],
                             collapse=params.collapse, build_tree=False)
    return test_data, test_set

def evaluate(prediction, methods, params, test_data=None, test_set = None):
    '''
    evaluate the prediction based on a test data set.
    parameters:
    prediction   --   an instance of predict_flu
    methods      --   a list of methods. each method contains a name, an extension, and a node set
    params       --   a name space with parameters including test regions and flutype
    '''
    flutype = params.flutype
    if 'toy'!=params.test:
        year=params.year
    
    if test_data is None:
        test_data, test_set = make_test_set(prediction,params)
        
    # calculate the various distances between predictions and the test set
    distances = {}
    distances_epi = {}
    distances['average'] = test_data.mean_distance_to_set(prediction.data.allele_frequencies)
    if flutype=='H3N2':
        distances_epi['average'] = test_data.aa_distance_to_set(prediction.data,positions = sorted(PF.HA1_antigenic_sites))
    print 'average\t', np.round(distances['average'], 4),
    if flutype=='H3N2':
        print np.round(distances_epi['average'], 4)

    # get the best possible prediction, either by loading the precomputed pickle
    # or by going through the prediction data and picking the best one post-hoc
    if os.path.isfile('../data/'+flutype+'_optimal_sequences_nucleotides.pickle'):
        with open('../data/'+flutype+'_optimal_sequences_nucleotides.pickle') as infile:
            optimal_seqs = pickle.load(infile)
        if test_set['start'].year in optimal_seqs:
            distances['minimal']= test_data.mean_distance_to_sequence(optimal_seqs[test_set['start'].year])
    if 'minimal' not in distances:
        distances['minimal'] = min([test_data.mean_distance_to_sequence(seq) for seq in prediction.data.aln])

    # if the flutype is H3N2, repeat for amino acid distances on epitopes defined for H3N2    
    if flutype=='H3N2' and os.path.isfile('../data/'+flutype+'_optimal_sequences_epitope.pickle'):
        with open('../data/'+flutype+'_optimal_sequences_epitope.pickle') as infile:
            optimal_seqs_epi = pickle.load(infile)
        if test_set['start'].year in optimal_seqs_epi:
            distances_epi['minimal']= test_data.aa_distance_to_sequence(optimal_seqs_epi[test_set['start'].year].seq.translate(),
                                                                        positions = sorted(PF.HA1_antigenic_sites))
    if flutype=='H3N2' and 'minimal' not in distances_epi:
        distances_epi['minimal'] = np.nan

    print 'minimal\t',np.round(distances['minimal']/distances['average'], 4),
    if flutype=='H3N2':
        print '\t',np.round(distances_epi['minimal']/distances_epi['average'], 4)

    # calculate the results from the different prediction methods
    for method in methods:
        mname = method[0]+method[1]
        distances[mname] = test_data.mean_distance_to_sequence(prediction.best_node(method[0], nodes = method[2]).seq)
        if flutype=='H3N2':
            distances_epi[mname] = test_data.aa_distance_to_sequence(prediction.best_node(method[0], 
                                         nodes = method[2]).aa_seq,positions = sorted(PF.HA1_antigenic_sites))
        print mname,'\t',np.round(distances[mname]/distances['average'], 4),
        if flutype=='H3N2':
            print '\t',np.round(distances_epi[mname]/distances_epi['average'], 4)
    
    # try to load L&L predictionss
    if os.path.isfile('../data/'+flutype+'_L_L_predictions.pickle'):
        laessig_prediction = get_LL(params.flutype)
        if test_set['start'].year in laessig_prediction:
            distances['L&L']= test_data.mean_distance_to_sequence(laessig_prediction[test_set['start'].year])
            distances_epi['L&L']= test_data.aa_distance_to_sequence(laessig_prediction[test_set['start'].year].seq.translate(),
                                                                    positions = sorted(PF.HA1_antigenic_sites))
            print 'L&L prediction\t',np.round(distances['L&L']/distances['average'],4),'\t',np.round(distances_epi['L&L']/distances_epi['average'],4)
        else:
            distances['L&L'] = np.nan
            distances_epi['L&L'] = np.nan
    else:
        distances['L&L'] = np.nan
    return distances, distances_epi, test_data


def make_combined_data(prediction, test_data, otherseqsnames=None, collapse = False):
    assert prediction.data.aln.get_alignment_length()==test_data.aln.get_alignment_length(),\
        "predict_and_test: prediction and test alignment have different length"
    # combined data set needed for plotting
    seqname_list = [seq.name for seq in prediction.data.aln] + [seq.name for seq in test_data.aln]
    if otherseqsnames:
        for seqname in otherseqsnames:
            if seqname not in seqname_list:
                seqname_list.append(seqname)

    combined_data = PF.flu_alignment(test_data.aln_file_name, test_data.outgroup, test_data.annotation,
                                     cds = test_data.cds, seq_names = seqname_list, collapse=collapse)
    
    return combined_data

