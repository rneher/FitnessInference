'''
run_flu_prediction.py

submit script that runs the prediction_statistics.py with different arguments
'''

import sys,os
import subprocess


nreps = 10
flutype = 'H3N2'
collapse = False

for ri in xrange(5):
    for boost in [0,0.5,1]:
        for dscale in [0.5]: #1.0, 2.0, 3.0, 5.0]:
            for D in [0.2, 0.5]:
                for sample_size in [50,100,200]:
                    for year in range(1995,2014):
                        for pred_region in ['"asia,north america"']:
                            arg_list = ['--flutype',flutype, '--year', year, '--sample_size', sample_size, '--nreps',nreps,
                                        '--pred', pred_region, '--test', pred_region, '--boost', boost, 
                                        '--eps_branch', 1e-5, '--dscale', dscale,'--diffusion',D]
                            if collapse: arg_list+=['--collapse']
                            cmd = 'qsub -cwd -l h_rt=0:59:59 -l h_vmem=8G ../src/prediction_statistics.py '\
                                +" ".join(map(str, arg_list))
                            print cmd
                            #os.system(cmd)

