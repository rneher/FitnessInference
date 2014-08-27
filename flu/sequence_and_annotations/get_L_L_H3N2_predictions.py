#############################
#get_L_L_H3N2_predictions.tex
#
# this script identifies the strain names we obtained from Luksza and Laessig
# with sequences in our mutliple sequence alignment
# 
################################
import sys
sys.path.append('../src')
import predict_flu as PF
from datetime import date
from Bio import Phylo,AlignIO,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip,os, argparse

#parse the command line arguments
aln_fname = '../data/H3N2_HA1_all_years_filtered.fasta.gz'

with open('../data/annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)

with gzip.open(aln_fname) as infile:
    total_alignment = AlignIO.read(infile,'fasta')
seq_look_up = {c.name.split('|')[0].lower():c for c in total_alignment}

L_L={}
with open('../data/H3N2_L_L_predicted_vaccine_strains.dat') as L_L_file:
    for line in L_L_file:
        entries = line.strip().split()
        year=int(entries[0][:-1])
        strain_name = '_'.join(entries[1].split('_')[:-3]).lower().replace('-', '_')
        if strain_name in seq_look_up:
            L_L[year-1] = seq_look_up[strain_name]
        else:
            print "strain for year", year-1, "not found"
with open('../data/H3N2_L_L_predictions.pickle', 'w') as outfile:
    pickle.dump(L_L, outfile)
