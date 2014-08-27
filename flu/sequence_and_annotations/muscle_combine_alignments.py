import glob
import subprocess
base_name  ='../data/gisaid_H3N2_all_years_human_full_date_???_aligned'
suffix ='.fasta'
flist = glob.glob(base_name  +suffix)
while len(flist)>1:
    flist.sort()
    suffix = 'a'+suffix
    for fi in range(0,len(flist)-1,2):
        cmd = ['muscle',  '-profile', '-in1', flist[fi], '-in2', flist[fi+1], '-out', base_name.replace('???', format(fi,'03d'))+suffix]
        subprocess.call(cmd)
        print ' '.join(cmd)
    if len(flist)%2:
        cmd = ['muscle',  '-profile', '-in1', flist[-1], '-in2',base_name.replace('???', format(fi,'03d'))+suffix , '-out', base_name.replace('???', format(len(flist)-1,'03d'))+suffix]
        subprocess.call(cmd)
        print ' '.join(cmd)
        
    flist= glob.glob(base_name+suffix)
    
