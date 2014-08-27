import sys,os, glob
import subprocess
import numpy as np


L=2000
sim_prefix = '/20140820_seqs'
dt=100
valdt = 1
jobcount = 0
D=0.5
year_list = range(5000, 24800, dt)
dir_list = glob.glob('../data_new/N_10000_L_2*sdt_1*')
for dir_name in dir_list:
    for sample_size in [200]:
        for dscale in [0.5, 1.0, 2.0, 3.0]:
            try:
                fname = dir_name+sim_prefix+'_prediction_results_'+str(sample_size)\
                        +'_d_'+str(dscale)+'_D_'+str(D)+'_dt_'+str(valdt*100)+'.dat'
                #print fname
                A = np.loadtxt(fname)
                years_done = list(A[:,0])
                #print years_done
                print dir_name, len(year_list)-len(years_done)
            except:
                print dir_name, 'nothing done'
                years_done = []
            for year in year_list:
                if year not in years_done:
                    args =  ['--base_name',dir_name+sim_prefix, '--sample_size', sample_size, '--dscale', dscale, '--diffusion',D,'--gen', year,'--valdt', valdt]
                    cmd = 'qsub -cwd -l h_rt=0:59:59 -l h_vmem=8G ../src/analyze_toy_data.py '\
                          + ' '.join(map(str, args))
                    print cmd, jobcount
                    os.system(cmd)
                    jobcount+=1
                if jobcount>2999:
                    exit()
