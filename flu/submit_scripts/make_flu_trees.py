###########################
# make_flu_trees.py
#
# submit script that loops over years and parameters and submits a job to 
# analyze a year of flu evolution , make a prediction and plot trees
###########################
import sys,os
import subprocess

flutype = 'H3N2'

for D in [0.5]:
    for boost in [0]:
        for dscale in [3.0,5.0]:
            for sample_size in [50,100,200]:
                for year in range(1995,2014):
                    for pred_region in ['"asia,north america"']:
                        arg_list = ['--year', year, '--sample_size', sample_size,'--flutype',flutype,
                                    '--pred', pred_region, '--test', pred_region, '--boost', boost, 
                                    '--eps_branch', 1e-5, '--dscale', dscale, '--diffusion', D, '--analysis', '--plot']
                        cmd = 'qsub -cwd -l h_rt=0:59:59 -l h_vmem=1G ../src/analyze_HA1.py '\
                            +" ".join(map(str, arg_list))
                        print cmd
                        os.system(cmd)

