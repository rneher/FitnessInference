'''
script that accepts the name of a fasta file as argument, reads all lines
and replaces 0 by A and 1 by T in the sequences
'''

import sys
import gzip

if len(sys.argv)==2:
    fname = sys.argv[1]
    foutname = '.'.join(fname.split('.')[:-1])+'_AT.fasta.gz'
    with gzip.open(foutname, 'w') as outfile:
        with open(fname, 'r') as infile:
            for line in infile:
                if line[0]=='>':
                    outfile.write(line)
                else:
                    outfile.write(line.replace('0','A').replace('1','T'))
