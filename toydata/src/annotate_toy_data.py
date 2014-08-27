import pickle, gzip,sys
from Bio import SeqIO

fitness_list = []
toy_data_annotations = {}
seq_file = sys.argv[1]
outgroup = [100000, None]
with gzip.open(seq_file, 'rb') as infile:
    for seq in SeqIO.parse(infile, 'fasta'):
        a={}
        a['colin']=False
        a['flubase']=False
        a['region']='toy'
        a['lat'] = 50
        a['lng'] = 0
        a['country'] = 'lego'
        a['date_info']='full_date'
        a['fitness'] = float(seq.name.split('_')[-1])
        fitness_list.append(a['fitness'])
        a['date'] = int(seq.name.split('_')[1])
        if a['date']<outgroup[0]:
            outgroup= [a['date'], seq]
        toy_data_annotations[seq.name]=a

with open('_'.join(seq_file.split('_')[:-1])+'_annotation.pickle', 'w') as outfile:
    pickle.dump(toy_data_annotations, outfile)

with open('_'.join(seq_file.split('_')[:-1])+'_outgroup.fa', 'w') as outfile:
    SeqIO.write(outgroup[1], outfile, 'fasta')

