'''
run_missing_flu_prediction.py

looks for aborted simulation (e.g. due to exceeded time limits) and resubmits them.

'''
import sys,os,time
import subprocess
import numpy as np
analysis_folder = 'analysis_may_feb/'

collapse = False
flutype = 'H3N2'
target_rep = 50
for boost in [0.0, 0.5, 1.0]:
    for dscale in [0.5]: #1.0, 2.0, 3.0, 5.0]:
        for D in [0.2, 0.5]:
            for sample_size in [50,100, 200]:
                for year in range(1995,2014):
                    for pred_region in ['"asia,north america"']:
                        base_name = flutype+'_'+str(year)+'_trees_asia_NA_'
                        name_mod = 'ssize_'+str(sample_size)+'_b_'+str(boost)+'_d_'+str(dscale)+'_D_'+str(D)
                        if collapse: name_mod+='_collapse'
                        fname_nuc = '_'.join(map(str,['../'+analysis_folder+flutype,'year',year, 'pred']+ 
                                                 ["-".join(pred_region[1:-1].split(','))]+['test']+
                                                 ['-'.join(pred_region[1:-1].split(',')),name_mod])).replace(' ', '')+'_nuc.dat'

                        try:
                            tmp_results = np.loadtxt(fname_nuc)
                            nreps = target_rep-tmp_results.shape[0]
                        except:
                            nreps=target_rep
                        print fname_nuc, nreps
                        if nreps>0:
                            print nreps
                            dn = 10
                            while nreps>0:
                                arg_list = ['--flutype',flutype, '--year', year, '--sample_size', sample_size, '--nreps',min(dn,nreps),
                                            '--pred', pred_region, '--test', pred_region, '--boost', boost, 
                                            '--eps_branch', 1e-5, '--dscale', dscale,'--diffusion',D]
                                if collapse: arg_list+=['--collapse']
                                cmd = 'qsub -cwd -l h_rt=0:59:59 -l h_vmem=1G ../src/prediction_statistics.py '\
                                      +" ".join(map(str, arg_list))
                                print cmd
                                os.system(cmd)
                                #time.sleep(.1)
                                nreps-=dn
