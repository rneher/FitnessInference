import sys,os,shutil, glob
import subprocess

prefix = '/20140418_'
L=10000
sigma = 0.1
bene_dele_ratio = 1
dir_list= glob.glob('../data_new/N_20000_L_10000*')
for dir_name in dir_list:
    if not os.path.isfile(dir_name+prefix+'seqs_annotation.pickle'):
        L,N,nflip,mu,sdt = map(float, dir_name.split('/')[-1].split('_')[1::2])
        if sdt == 100:
            ssize = 200
        elif sdt ==1:
            ssize = 2
        arg_list = [L,N, ssize, sdt, sigma, nflip, mu, bene_dele_ratio]
        cmd = 'qsub -cwd -l h_rt=10:59:59 -l h_vmem=8G submit_script.py '\
              +" ".join(map(str, arg_list))
        print cmd
        os.system(cmd)
    else:
        print dir_name,'OK'
        
