##################
#find_duplicates.py
#script to filter sequences such that a single place and date has only unique sequences.
##################
from datetime import date
from Bio import Phylo,SeqIO
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import glob,pickle,gzip
from collections import defaultdict

with open('annotations.pickle', 'r') as infile:
    annotation = pickle.load(infile)

seqs_1968 = []
seqs = defaultdict(list)
with gzip.open('H3N2_HA1_all_years_filtered.fasta.gz', 'rb') as infile:
    for seq in SeqIO.parse(infile, 'fasta'):
        if annotation[seq.name]['date_info']=='full_date':
            strain = seq.name.split('|')[0]
            place = strain.split('/')[1]
            seqs[(place, annotation[seq.name]['date'], str(seq.seq))].append(seq)
        elif 'Hong_Kong/68' in seq.name:
            seqs_1968.append(seq)
with gzip.open('H3N2_HA1_all_years_filtered_unique.fasta.gz', 'wb') as infile:
    for seq in seqs_1968:
        SeqIO.write(seq, infile, 'fasta')
    for k, seq_list in seqs.iteritems():
        SeqIO.write(seq_list[0], infile, 'fasta')

