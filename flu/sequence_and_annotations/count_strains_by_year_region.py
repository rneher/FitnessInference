###################
# count_strains_by_year_region.py
#script that prints a table of sequence counts in different years and geographic
#regions.
###################
from StringIO import StringIO
from socket import gethostname
import sys
sys.path.append('../../prediction_src')
import predict_flu as flu
from collections import defaultdict
import gzip,pickle
from Bio import SeqIO

flu_type = 'H3N2_gisaid'
flu_db=False
if flu_db:
    with open('../data/annotations.pickle', 'r') as infile:
        annotations = pickle.load(infile)

    seqs_by_year = defaultdict(list)
    flutype = 'H3N2'
    with gzip.open('../data/'+flutype+'_HA1_all_years_filtered.fasta.gz', 'r') as infile:
        for seq in SeqIO.parse(infile, 'fasta'):
            seq_name =  seq.name #.split('|')[0]
            if annotations[seq_name]['date_info'] in ['full_date']:
                seqs_by_year[(annotations[seq_name]['date'].year, 
                              flu.coordinates_to_region(annotations[seq_name]['lat'],
                                                        annotations[seq_name]['lng'])
                              )].append(seq)
else:
    with open('../data/annotations_gisaid.pickle', 'r') as infile:
        annotations = pickle.load(infile)

    seqs_by_year = defaultdict(list)
    with gzip.open('../data/gisaid_H3N2_all_years_human_full_date_aligned_trimmed.fasta.gz', 'r') as infile:
        for seq in SeqIO.parse(infile, 'fasta'):
            seq_name =  seq.name #.split('|')[
            if seq_name in annotations and annotations[seq_name]['date_info'] in ['full_date']:
                seqs_by_year[(annotations[seq_name]['date'].year, 
                              flu.coordinates_to_region(annotations[seq_name]['lat'],
                                                        annotations[seq_name]['lng'])
                              )].append(seq)



counts = {}
regions = ['north america', 'south america', 'europe', 'asia', 'oceania']
for year in range(1970, 2015):
    counts[year] = {}
    for region in regions:
        label = (year, region)
        if label in seqs_by_year:
            counts[year][region]=len(seqs_by_year[label])
        else:
            counts[year][region]=0

with open('../data/'+flu_type+'_seqs_by_year_and_region.txt', 'w') as outfile:
    print 'year',
    outfile.write('year')
    for region in regions:
        print '\t',region,
        outfile.write('\t'+region)
    print
    outfile.write('\n')

    seq_sum = 0
    for year in sorted(counts.keys()):
        print year,
        outfile.write(str(year))
        for region in regions:
            print '\t', counts[year][region],
            outfile.write('\t'+str(counts[year][region]))
            seq_sum+=counts[year][region]
        print
        outfile.write('\n')
