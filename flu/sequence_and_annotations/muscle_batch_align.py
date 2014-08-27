import glob
import os

flist = glob.glob('../data/gisaid_H3N2_all_years_human_full_date_???.fasta')

for fname in flist:
    cmd = ['qsub', '-cwd -l h_vmem=8G -l h_rt=00:59:59', 'muscle_script.sh',  fname, fname[:-6]+'_aligned.fasta']
    os.system(' '.join(cmd))
    print ' '.join(cmd)
