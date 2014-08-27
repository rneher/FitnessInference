from Bio import SeqIO
from datetime import  date

fname = '../data/gisaid_H3N2_all_years_human.fasta'
all_seqs = []

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

with open(fname, 'r') as seqfile:
    for seq in SeqIO.parse(seqfile, 'fasta'):
        collection_date = parse_gisaid_date(seq.description.split('|')[1].strip())
        if collection_date is not None:
            all_seqs.append([seq, collection_date])

all_seqs.sort(key = lambda x:x[1])
outfname = '../data/gisaid_H3N2_all_years_human_full_date_'
step = 500
for si in range(len(all_seqs)/step+1):
    with open(outfname+format(si,'03d')+'.fasta', 'w') as outfile:
        for seq, d in all_seqs[si*step:(si+1)*step]:
            SeqIO.write(seq, outfile, 'fasta')

                              

