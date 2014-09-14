import glob,sys
import numpy as np
sys.path.append('../../flu/src')
import test_flu_prediction as test_flu
import matplotlib.pyplot as plt
import analysis_utils_toy_data as AU

file_formats = ['.svg', '.pdf']
plt.rcParams.update(test_flu.mpl_params)


line_styles = ['-', '--', '-.']
cols = ['b', 'r', 'g', 'c', 'm', 'k', 'y']
cols+=cols
figure_folder = '../figures/'
data_dir = '../data_new'
prefix= '/20140820_'

N_list = [20000] #,20000]
mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6]
#nflip_list = [0.02,0.04, 0.08, 0.16]
gamma_list = [1.0] #, 2.0,3.0, 5.0]
omega_list = [0.3]
nflip_list = [0.04, 0.08]
sdt_list = [1,100]  #determines whether 2 genomes are sampled every generation, or 200 every 100 gen
pred, norm_pred, run_stats = AU.load_prediction_data(prefix, N_list, mu_list, nflip_list,
                                                     sdt_list, return_mean=True)

valdt = 200
ssize = 200
D = 0.2
L=2000
mean_fitness_true_fitness_spearman_i = -4

for gamma in gamma_list:
    for omega in omega_list:
        for sdt in [1,100]:
            plt.figure(figsize= (10,6))
            ax = plt.subplot(111)
            #plt.title(r'\omega='+str(omega)+',\;dt='+str(sdt)+'$')    
            for di,D in enumerate([0.2, 0.5]):
                pred_label = ssize, gamma, D, omega, valdt
        
                ### PLOT FITNESS CORRELATION VS PAIRWISE DIVERSITY ###
                for ni, N in enumerate(N_list):
                    for fi, nflip in enumerate(nflip_list):
                        plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                                 [pred[(N,mu,nflip,sdt)+pred_label][0][mean_fitness_true_fitness_spearman_i] for mu in mu_list],
                                 [pred[(N,mu,nflip,sdt)+pred_label][2][mean_fitness_true_fitness_spearman_i] for mu in mu_list],
                                 c=cols[fi], ls=line_styles[di], label = '$n_A = '+str(nflip)+',\;\Gamma='+str(D)+'$', lw=2)
                    
            plt.ylabel("Spearman's correlation")
            plt.xlabel('average pairwise distance')
            plt.xscale('log')
            plt.legend(loc=4)
            #add panel label 
            plt.text(0.02,0.9,'Fig.~2-S1', transform = plt.gca().transAxes, fontsize = 20)
            plt.xlim([0.5, 200])
            for ff in file_formats:
                plt.savefig(figure_folder+'Fig2_S1_pairwise_diversity_vs_predictability_sdt_'+str(sdt)+'_gamma_'
                        +str(gamma)+'_valdt_'+str(valdt)+ff)
    
    
            ### PLOT prediction_success VS PAIRWISE DIVERSITY ###
            plt.figure(figsize= (10,6))
            ax = plt.subplot(111)
            #plt.title(r'$\gamma='+str(gamma)+',\;\omega='+str(omega)+',\;dt='+str(sdt)+'$')    
    
            for ni, N in enumerate(N_list):
                for fi, nflip in enumerate(nflip_list):
                    for di,D in enumerate([0.2, 0.5]):
                        plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                                 [norm_pred[(N,mu,nflip,sdt)+pred_label][0][1] for mu in mu_list],
                                 [norm_pred[(N,mu,nflip,sdt)+pred_label][2][1] for mu in mu_list],
                                 c=cols[fi], ls=line_styles[di], label = '$n_A = '+str(nflip)+'$')
    
            plt.ylabel(r'distance $\bar{d}$ to future populations')
            plt.xlabel('average pairwise distance')
            #add panel label 
            plt.text(0.02,0.9,'Fig.~2-S2', transform = plt.gca().transAxes, fontsize = 20)
            plt.xscale('log')
            plt.legend(loc=1)
            plt.xlim([0.5, 200])
    
            for ff in file_formats:
                plt.savefig(figure_folder+'Fig2_S2_pairwise_diversity_vs_distance_sdt_'+str(sdt)+'_gamma_'+str(gamma)+'_valdt_'+str(valdt)+ff)
            #plt.close()


## plot gamma versus the number of predictions that are worse than random
# reload the data without averaging over the different realizations.
pred, norm_pred, run_stats = AU.load_prediction_data(prefix, N_list, mu_list, nflip_list,
                                                     sdt_list, return_mean=False)


#for sdt in [100]:
if len(gamma_list)>1:
    for omega in omega_list:
    
            ### PLOT FITNESS CORRELATION VS DSCALE ###
        plt.figure(figsize= (10,6))
        ax = plt.subplot(111)
        #plt.title(r'$\omega='+str(omega)+',\;dt='+str(sdt)+'$')    
        for mi,mu in enumerate(mu_list[3:6]):
            for ni, N in enumerate(N_list):
                for fi, nflip in enumerate(nflip_list[:]):
                    if mi==0:
                        label_str = r'$n_A ='+str(nflip)+'$'
                    else:
                        label_str = None
                    plt.plot(gamma_list, [np.mean(pred[(N,mu,nflip,sdt)+(ssize, gamma, D, omega, valdt)][:,0]<
                                                  pred[(N,mu,nflip,sdt)+(ssize, gamma, D, omega, valdt)][:,3])
                                          for gamma in gamma_list], lw=2, marker='o', markersize=10,
                                    ls=line_styles[mi], c=cols[fi], label = label_str)

        #plt.xscale('log')
        #add panel label 
        plt.text(0.02,0.9,'Fig.~2-S3', transform = plt.gca().transAxes, fontsize = 20)
        plt.xlim([0.0, 5.5])
        plt.ylabel('worse than random (out of 100)')
        plt.xlabel(r'time rescaling $\gamma$')
        plt.legend(loc=1,numpoints=1)
        for ff in file_formats:
            plt.savefig(figure_folder+'Fig2_S3_gamma_vs_predictability_sdt_'+str(sdt)+'_D_'+str(D)+'_w_'+str(omega)+'_valdt_'+str(valdt)+ff)

