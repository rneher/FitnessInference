
import sys,os,shutil
import subprocess

L=2000
sigma = 0.03
bene_dele_ratio = 1

for N in [10000, 20000]:
    for sdt,ssize in [(1,2)]:
        for mu in [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6]:
            for nflip in [0.01,0.02,0.04, 0.08, 0.16]:
                try:
                    dirname = '_'.join(map(str,['../data_new/N', N, 'L', L, 'nflip',nflip,'mu',mu,'sdt', sdt]))
                    os.mkdir(dirname)
                    os.mkdir(dirname+'/src')
                    shutil.copy('../src/flusim.cpp', dirname+'/src')
                    shutil.copy('../src/annotate_toy_data.py', dirname+'/src')
                except OSError, e:
                    print "cannot make directory", e
                arg_list = [L,N, ssize, sdt, sigma, nflip, mu, bene_dele_ratio]
                cmd = 'qsub -cwd -l h_rt=10:59:59 -l h_vmem=8G submit_script.py '\
                    +" ".join(map(str, arg_list))
                print cmd
                os.system(cmd)

