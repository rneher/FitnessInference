'''
makes a dictionary of geographic location found in flubase name strings.
these name strings are submitted to google geocoder to obtain full
geographic information including longitude and latidue. based on this,
regions like south america are defined
'''


from Bio import SeqIO
import pickle
import time,gzip

fname = '../data/H3N2_HA1_2012_2014_flubase_filtered.fasta.gz'
#fname = '../data/H1N1_HA1_all_years_filtered.fasta.gz'
#fname = '../data/gisaid_H3N2_all_years_human_full_date_aligned_trimmed.fasta.gz'
places = set()

names = []
with gzip.open(fname, 'r') as infile:
    for seq in SeqIO.parse(infile, 'fasta'):
        names.append(seq.name)
        #country = seq.name.split('|')[3]
        country = ''
        places.add((seq.name.split('/')[1],country))

# append to an existing pickle -- this is useful since geocoders have a quota
# one global file for all flu strains should be enough
try:
    with open('../data/flubase_places.pickle', 'r') as infile:
        place_to_coordinates = pickle.load(infile)
except:
    place_to_coordinates = {}

from geopy import geocoders
g = geocoders.GoogleV3()
g.timeout=10
for place, country in places:
    if place not in place_to_coordinates:
        loc = g.geocode(place.replace('_', ' ')+', '+country.replace('_', ' '))
        time.sleep(0.2)
        try:
            if loc:
                print place, country,loc
                country = loc[0].split(',')[-1].strip()
                country = country.encode('ascii', 'replace')
                location = loc[1]
            else:
                print place, country
                location = ('nan', 'nan')
                country = 'unknown'

            place_to_coordinates[place] = {'country':country.lower(), 
                                           'lat':location[0], 'lng':location[1]}
        except:
            print "ERROR",place

# save pickle back to the same file with new places added
try:
    with open('../data/flubase_places.pickle', 'w') as outfile:
        pickle.dump(place_to_coordinates,outfile)
except:
    print "can't save places_pickle"
