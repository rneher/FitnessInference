#!/ebio/ag-neher/share/programs/bin/python
#!/usr/bin/python

import sys,os,glob

# run the compiled C++ program to generate the data
executable = '/ebio/ag-neher/share/users/rneher/FluPrediction_code/toydata/src/flu'
cmd = executable+' '+' '.join(sys.argv[1:])
print cmd
os.system(cmd)

# reconstruct the destination directory and translate the output to an AT format rather than binary 01010
outdirname = '_'.join(map(str,['../data_new/N', sys.argv[2], 'L', sys.argv[1],
                               'nflip',sys.argv[6],'mu',sys.argv[-2],'sdt', sys.argv[4]]))
flist = glob.glob(outdirname+'/*seqs.fa')
print outdirname ,flist
for fname in flist:
    cmd = 'python ../src/translate_binary_fasta.py '+fname
    os.system(cmd)
    cmd = 'rm '+fname
    os.system(cmd)

# make an annotation file.
flist = glob.glob(outdirname+'/*_AT.fasta.gz')
print outdirname ,flist
for fname in flist:
    cmd = 'python ../src/annotate_toy_data.py '+fname
    os.system(cmd)
