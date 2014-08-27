'''
remove_outliers_from_seqs.py

open a file with a list of strain names deemed outliers, opens the file
with all sequences, and prints those to file that are not in the outliers file.
'''
from Bio import Phylo, AlignIO, SeqIO

HA_seqs = AlignIO.read('../data/H1N1_HA1_all_years_trimmed.fasta', 'fasta')

outliers = []
with open('../data/H1N1_HA1_all_years_outliers.txt', 'r') as infile:
    for line in infile:
        outliers.append(line.strip().split('|')[0])


with open('../data/H1N1_HA1_all_years_filtered.fasta', 'w') as outfile:
    for seq in HA_seqs:
        if seq.name.split('|')[0] not in outliers:
            seq.name = '/'.join(seq.name.split('/')[:-1])
            SeqIO.write(seq, outfile, 'fasta')

