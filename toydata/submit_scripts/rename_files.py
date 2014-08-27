import sys,os
import subprocess
import time, glob

L=2000
sim_prefix = '/20140820_seqs'
dt=200
jobcount = 0
# inner loops = 2*2*4*1*3*4 = 2**7 * 3 approx 400
n=20
i=0
D=0.5
istep=0
first_gen = 5000
last_gen = 25000 #min(5000 + istep*(i+1)*dt, 24800)
dir_list= glob.glob('../data_new/N_10000_L_2*_sdt_1*')
for year in xrange(first_gen, last_gen, n*dt):
    for dir_name in dir_list:
        flist = glob.glob(dir_name+sim_prefix+'*prediction_results*D_?.?.dat')
        for fname in flist:
            cmd = 'mv '+  fname +' ' + fname[:-4]+'_dt_100.dat'
            os.system(cmd)
            print cmd
              
