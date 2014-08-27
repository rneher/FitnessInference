##############################
#get_optimal_seqs.py
#
# this script identifies the best possible sequence for each year and writes them
# into a pickle file for easy reuse by other scripts
#################################
import sys
sys.path.append('../../prediction_src')
import predict_flu as PF
from datetime import date
from Bio import Phylo,AlignIO,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os, argparse

sample_size = 10000
years = range(2013,2014)
prediction_results = []
cds = {'begin':0, 'end':987, 'pad':0}

#parse the command line arguments
parser = argparse.ArgumentParser(description="predict the flu")
parser.add_argument('--flutype', type = str, default = 'H3N2', help = 'influenza strain')
parser.add_argument('--pred', default='asia,north america', type=str, help='regions to predict')
parser.add_argument('--test', default='asia,north america', type=str, help='regions to test')
params=parser.parse_args()
prediction_regions = params.pred.split(',')
test_regions = params.test.split(',')
flutype = params.flutype
aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'


with open('../data/'+flutype+'_annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)
outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')

prediction_set={}
test_set = {}
if os.path.isfile('../data/'+flutype+'_optimal_sequences_nucleotides.pickle'):
    with open('../data/'+flutype+'_optimal_sequences_nucleotides.pickle', 'r') as infile:
        optimal_sequences_nucs = pickle.load(infile)
else:
    optimal_sequences_nucs = {}
    
if os.path.isfile('../data/'+flutype+'_optimal_sequences_epitopes.pickle'):
    with open('../data/'+flutype+'_optimal_sequences_epitope.pickle', 'r') as infile:
        optimal_sequences_epi = pickle.load(infile)
else:
    optimal_sequences_epi = {}

for year in years:
    prediction_set={'start': date(year-2, 4,1), 'stop':date(year, 2,28),
                    'regions':prediction_regions, 'sample_size':sample_size}
    test_set = {'start':date(year, 10,1), 'stop':date(year+1, 3,31),
                'regions':test_regions, 'sample_size':50}

    prediction_data = PF.flu_data_set(aln_fname,outgroup,annotation,cds)
    prediction_data.select_subset([[prediction_set['start'], prediction_set['stop'],
                                    [reg],prediction_set['sample_size']] for reg in prediction_set['regions']])

    test_data = PF.flu_data_set(aln_fname,outgroup,annotation,cds)
    test_data.select_subset([[test_set['start'], test_set['stop'],
                            [reg],test_set['sample_size']] for reg in test_set['regions']])

    optimal_sequences_nucs[year] = min([(test_data.mean_distance_to_sequence(seq),seq) for seq in prediction_data.aln])[1]
    optimal_sequences_epi[year] = min([(test_data.aa_distance_to_sequence(seq.seq.translate(),positions = PF.HA1_antigenic_sites),seq)
                                        for seq in prediction_data.aln if '-' not in seq.seq])[1]


with open('../data/'+flutype+'_optimal_sequences_nucleotides.pickle', 'w') as outfile:
    pickle.dump(optimal_sequences_nucs, outfile)

with open('../data/'+flutype+'_optimal_sequences_epitope.pickle', 'w') as outfile:
    pickle.dump(optimal_sequences_epi, outfile)
