import sys,os
import subprocess
import time, glob

L=2000
sim_prefix = '/20140820_seqs'
dt=100
valdt = 1
jobcount = 0
# inner loops = 2*2*4*1*3*4 = 2**7 * 3 approx 400
n=10
i=0
D=0.5
istep=0
first_gen = 5000
last_gen = 25000 #min(5000 + istep*(i+1)*dt, 24800)
dir_list= glob.glob('../data_new/N_10000_L_2*_sdt_1*')
for year in xrange(first_gen, last_gen, n*dt):
    for dir_name in dir_list:
        #cmd = 'rm '+dir_name+sim_prefix+'*prediction_results*dt_?00.dat'
        #os.system(cmd)
        if True:
            for dscale in [0.5]: #, 1.0,2.0,3.0]:
                for sample_size in [200]:
                    args = ['--base_name',dir_name+sim_prefix, '--sample_size', sample_size, '--dscale', dscale, '--diffusion',D,'--gen', year '--valdt', valdt, n, dt]
                    cmd = 'qsub -cwd -l h_rt=0:59:59 -l h_vmem=8G analyze_multiple_toy_data.py '\
                          + ' '.join(map(str, args))

                    print cmd
                    os.system(cmd)
                    jobcount+=1
                    if jobcount>3000:
                        exit()
