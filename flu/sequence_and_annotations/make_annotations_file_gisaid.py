'''
loops over a sequence aligment and attempts to parse geographic location and
date from the name. this information is stored as a dictionary with sequence name
as key and saved in a file annotations.pickle. this dictionary can be used cumulatively
for multiple sequence input files.
'''
from Bio import Phylo, AlignIO, SeqIO
from datetime import date
import pickle, gzip

flu_type = 'H3N2_gisaid'
fname = '../data/'+flu_type+'_HA1_all_years_filtered.fasta.gz'

# open place list. this is data set unspecific
with open('../data/flubase_places.pickle', 'r') as infile:
    places_to_coordinates = pickle.load(infile)

# read alignment
with gzip.open(fname) as infile:
    HA_seqs = AlignIO.read(infile, 'fasta')

HA_seqs_dict = {leaf.name.split('|')[0]:leaf for leaf in HA_seqs}
full_date, year_only, no_date_info = 'full_date', 'year_only', 'no_date_info'

def parse_gisaid_date(date_str):
    if len(date_str.split('-'))==3:
        year, month, day = map(int, date_str.split('-'))
        return date(year =year, month=month, day=day)
    elif len(date_str.split('-'))==2:
        year, month = map(int, date_str.split()[0].split('-'))
        print date_str
        return date(year =year, month=month, day=15)
    else:
        None

# add to existing file of files exists. we only need one file for all flu data sets. 
try:
    with open('../data/'+flu_type+'_annotations.pickle', 'r') as infile:
w        annotations = pickle.load(infile)
except:
    annotations = {}
    print "made new annotations file"
for seq in HA_seqs:
    try:
        
        if seq.name not in annotations:
            if 'gisaid' in flu_type:
                seq_date = parse_gisaid_date(seq.description.split('|')[1].strip())
            else:
                seq_date = parse_gisaid_date(seq.description.split('|')[1].strip())
            if seq_date is not None:
                date_info = 'full_date'
                place = seq.name.split('/')[1]
                if place in places_to_coordinates:
                    annotations[seq.name]={}
                    annotations[seq.name]['colin']=False
                    annotations[seq.name]['flubase']=False
                    annotations[seq.name]['gisaid']=True            
                    annotations[seq.name]['date'] = seq_date
                    annotations[seq.name]['date_info']=date_info
                    annotations[seq.name]['country'] = places_to_coordinates[place]['country'].upper()
                    annotations[seq.name]['lat'] =places_to_coordinates[place]['lat']
                    annotations[seq.name]['lng'] =places_to_coordinates[place]['lng']
                    annotations[seq.name]['coordinates'] = (places_to_coordinates[place]['lat'],
                                                        places_to_coordinates[place]['lng'])
                else:
                    print "no place for:", seq.name
    except:
        print "no annotation for", seq.name

with open('../data/'+flu_type+'_annotations.pickle', 'w') as outfile:
    pickle.dump(annotations, outfile)


