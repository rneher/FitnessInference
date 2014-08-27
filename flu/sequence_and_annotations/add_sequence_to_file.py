################
#script that filters sequences that contain ambiguous nucleotides or deletions
#
################
from Bio import Phylo, AlignIO, SeqIO

file1 = '../data/H3N2_HA1_all_years_filtered.fasta'
file2 = '../data/H3N2_HA1_2012_2014_flubase_filtered.fasta'

HA_seqs_file1 = AlignIO.read(file1,'fasta')
HA_seqs_file2 = AlignIO.read(file2,'fasta')

seq_names = []
for leaf in HA_seqs_file1:
    seq_names.append('|'.join(leaf.name.split('|')[:1]))


with open(file1, 'a') as outfile:
    count_present = 0
    count_added = 0
    for seq in HA_seqs_file2:
        s_name = '|'.join(seq.name.split('|')[:1])
        if s_name in seq_names:
            print 'PRESENT',seq.name
            count_present+=1
        else:
            print 'ADDED', seq.name
            #SeqIO.write(seq, outfile, 'fasta')
            count_added+=1
print 'present', count_present
print 'added', count_added

        
