################
#script that filters sequences that contain ambiguous nucleotides or deletions
#
################
from Bio import Phylo, AlignIO, SeqIO

prefix = 'H3N2_HA1_2012_2014_flubase'
#HA_seqs = AlignIO.read('../data/H1N1_HA1_all_years_filtered.fasta', 'fasta')
HA_seqs = AlignIO.read('../data/'+prefix+'_trimmed.fasta', 'fasta')
for leaf in HA_seqs:
    leaf.name = '/'.join(leaf.name.split('/')[:-1])

low_q_seqs = []
good_seqs = []
for seq in HA_seqs:
    if sum([seq.seq.count(c) for c in 'HKMNRVWSXYKF'])<4 and seq.seq.count('-')<4:
        good_seqs.append(seq)
    else:
        low_q_seqs.append(seq)

with open('../data/'+prefix+'_low_quality_strains.txt', 'w') as outfile:
    for seq in low_q_seqs:
        outfile.write(seq.name+'\n')

with open('../data/'+prefix+'_filtered.fasta', 'w') as outfile:
    for seq in good_seqs:
        SeqIO.write(seq, outfile, 'fasta')
        
