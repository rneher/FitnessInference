'''
loops over a sequence aligment and attempts to parse geographic location and
date from the name. this information is stored as a dictionary with sequence name
as key and saved in a file annotations.pickle. this dictionary can be used cumulatively
for multiple sequence input files.
'''
from Bio import Phylo, AlignIO, SeqIO
from datetime import date
import pickle, gzip

fname = '../data/H3N2_HA1_2012_2014_flubase_filtered.fasta.gz'
#fname = '../data/H1N1_HA1_all_years_filtered.fasta.gz'


# open place list. this is data set unspecific
with open('../data/flubase_places.pickle', 'r') as infile:
    places_to_coordinates = pickle.load(infile)

# read alignment
with gzip.open(fname) as infile:
    HA_seqs = AlignIO.read(infile, 'fasta')

HA_seqs_dict = {leaf.name.split('|')[0]:leaf for leaf in HA_seqs}
full_date, year_only, no_date_info = 'full_date', 'year_only', 'no_date_info'

def parse_date_string(dstr):
    '''
    function that takes a flubase date and attempts conversion to a python date
    structure. checks for mm/dd/yyyy
    if xx/xx structure is found, it assumed mm/yy
    if total length is 4 without /, assumes yyyy
    if total length is 2 without /, assumes yy
    '''
    try:
        if len(dstr.split('/'))==3:
            entries = map(int, dstr.split('/'))
            return date(year=entries[2], month = entries[0], day = entries[1]), full_date
        elif len(dstr.split('/'))==2:
            entries = map(int, dstr.split('/'))
            return date(year=entries[1], month = entries[0], day = 15), full_date
        elif len(dstr)==4:
            return date(year=int(dstr), month = 7, day = 1), year_only
        elif len(dstr)==2:
            year = int(dstr)
            if year<20:
                year+=2000
            else:
                year+=1900
            return date(year=year, month = 7, day = 1), year_only
        else:
            return None, no_date_info
    except:
        return None, no_date_info


# add to existing file of files exists. we only need one file for all flu data sets. 
try:
    with open('../data/H3N2_annotations.pickle', 'r') as infile:
        annotations = pickle.load(infile)
except:
    annotations = {}
    print "made new annotations file"
for seq in HA_seqs:
    try:
        if seq.name not in annotations:
            annotations[seq.name]={}
            annotations[seq.name]['colin']=False
            annotations[seq.name]['flubase']=True
            seq_date, date_info = parse_date_string(seq.name.split('|')[2])
            annotations[seq.name]['date'] = seq_date
            annotations[seq.name]['date_info']=date_info
            place = seq.name.split('/')[1]
            country = seq.name.split('|')[3]
            if place in places_to_coordinates:
                annotations[seq.name]['country'] = places_to_coordinates[place]['country'].upper()
                annotations[seq.name]['lat'] =places_to_coordinates[place]['lat']
                annotations[seq.name]['lng'] =places_to_coordinates[place]['lng']
                annotations[seq.name]['coordinates'] = (places_to_coordinates[place]['lat'],
                                                    places_to_coordinates[place]['lng'])
            else:
                print "no place for:", seq.name
    except:
        print "no annotation for", seq.name

with open('../data/H3N2_annotations.pickle', 'w') as outfile:
    pickle.dump(annotations, outfile)


