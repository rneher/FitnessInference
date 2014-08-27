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


import pickle,os
import predict_flu as PF
import numpy as np

def predict(aln_fname, outgroup, annotation, methods, prediction_set, cds,
                     time_bins=None, subsample_factor = 1.0, boost = 0.0,
                     eps_branch_length = 1e-5, collapse = False, dscale = 3.0, D=0.5, pseudo_count=5):
    '''
    simple wrapper for the predict_flu class. 
    produces a flu data set given the criteria and predicts using this
    data set and the parameters passed
    '''
    # load data for prediction and build tree
    prediction_data = PF.flu_alignment(aln_fname,outgroup,annotation,cds=cds, subsample_factor = subsample_factor,
                                      criteria = [[prediction_set['start'], prediction_set['stop'],
                                    [reg],prediction_set['sample_size']] for reg in prediction_set['regions']], collapse=collapse)

    # PREDICT
    prediction = PF.flu_ranking(prediction_data, eps_branch_length = eps_branch_length, boost = boost,
                                time_bins = time_bins, methods=methods, D=D,
                                 distance_scale = dscale,pseudo_count = pseudo_count)
    prediction.predict()
    print "prediction done"
    return prediction


def evaluate(prediction, methods, test_set, aln_fname, outgroup, annotation, cds,flutype='H3N2', collapse=False):
    '''
    evaluate the prediction based on a test data set.
    parameters:
    prediction   --   an instance of predict_flu
    methods      --   a list of methods. each method contains a name, an extension, and a node set
    test_set     --   set of start,stop data, regions and sample size to test on
    aln_fname    --   alignment containing the test_data
    outgroup     --   outgroup to root the tree. 
    cds          --   coding sequence coordinates {begin:??, end:??, pad:??}
    '''
    # load data to test
    test_data = PF.flu_alignment(aln_fname,outgroup,annotation,cds=cds,
                             criteria = [[test_set['start'], test_set['stop'],
                            [reg],test_set['sample_size']] for reg in test_set['regions']], collapse=collapse)
    
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
        with open('../data/'+flutype+'_L_L_predictions.pickle') as infile:
            laessig_prediction = pickle.load(infile)
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

