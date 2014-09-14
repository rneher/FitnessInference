import glob,sys
import numpy as np
sys.path.append('../../flu/src')
import test_flu_prediction as test_flu
import matplotlib.pyplot as plt
from matplotlib import cm
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
mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6] #,128e-6]
#nflip_list = [0.02,0.04, 0.08, 0.16]
offset = 2
# multiply by 2 to transform to pairwise distance
m_list = 2.0**np.arange(-6,4, 1) * 0.5
m_to_plot = m_list[offset:-2]
nflip = 0.08 #, 0.16]
sdt_list = [1,100]  #determines whether 2 genomes are sampled every generation, or 200 every 100 gen
pred, norm_pred, run_stats, corrcoef = AU.load_prediction_data(prefix, N_list, mu_list, [nflip],
                                                     sdt_list, return_mean=True, polarizer=True)

pred_I, norm_pred_I, run_stats_I = AU.load_prediction_data(prefix, N_list, mu_list, [nflip],
                                                     sdt_list, return_mean=True, polarizer=False)

D,gamma,omega = 0.2,1.0,0.3

valdt = 200
ssize = 200
L=2000
mean_fitness_true_fitness_spearman_i = -4

###################################################################
### correlation vs pairwise diversity
###################################################################
for sdt in [1,100]:
    pred_label = (valdt,)
    pred_label_I = ssize, gamma, D, omega, valdt

    ### PLOT FITNESS CORRELATION VS PAIRWISE DIVERSITY ###
    plt.figure(figsize= (8,6))
    #plt.title('sdt = '+str(sdt))
    ax = plt.subplot(111)
    for ni, N in enumerate(N_list):
        for mi,m in enumerate(m_to_plot):
            print [corrcoef[(N,mu,nflip,sdt)+pred_label][0][mi+offset] for mu in mu_list]
#                if mi: label_str = None
#                else: label_str='$n_A = '+str(nflip)+'$'
            label_str=r'$\tau = '+str(m)+'$'
            plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                     [corrcoef[(N,mu,nflip,sdt)+pred_label][0][mi+offset] for mu in mu_list],
                     [corrcoef[(N,mu,nflip,sdt)+pred_label][2][mi+offset] for mu in mu_list],
                     c=cm.jet((mi+1.0)/(1+len(m_to_plot))), ls='-', label = label_str , lw=2)

        plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                 [pred_I[(N,mu,nflip,sdt)+pred_label_I][0][mean_fitness_true_fitness_spearman_i] for mu in mu_list],
                 [pred_I[(N,mu,nflip,sdt)+pred_label_I][2][mean_fitness_true_fitness_spearman_i] for mu in mu_list],
                 c='k', ls='--', label = 'Fitness inference', lw=2)
            
    plt.ylabel("Spearman's correlation")
    plt.xlabel('average pairwise distance')
    plt.xscale('log')
    plt.legend(loc=4)
    #add panel label 
    #plt.text(0.02,0.9,'Fig.~1-S1', transform = plt.gca().transAxes, fontsize = 20)
    plt.xlim([1.0, 100])
    plt.tight_layout()
    for ff in file_formats:
        plt.savefig(figure_folder+'Fig3_pairwise_diversity_vs_predictability_polarizer_sdt_'+str(sdt)+'_nflip_'+str(nflip)+'_valdt_'+str(valdt)+ff)


    ### PLOT prediction_success VS PAIRWISE DIVERSITY ###
    plt.figure(figsize= (8,6))
    ax = plt.subplot(111)
    #plt.title(r'$dt='+str(sdt)+'$')    

    for ni, N in enumerate(N_list):
        for mi,m  in enumerate(m_to_plot):
#                if mi: label_str = None
#                else: label_str='$n_A = '+str(nflip)+'$'
            label_str=r'$\tau = '+str(m)+'$'
            plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                     [norm_pred[(N,mu,nflip,sdt)+pred_label][0][2+mi] for mu in mu_list],
                     [norm_pred[(N,mu,nflip,sdt)+pred_label][2][2+mi] for mu in mu_list],
                     c=cm.jet((mi+1.0)/(1+len(m_to_plot))), ls='-',  label =label_str ,lw=2)

        plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                 [norm_pred_I[(N,mu,nflip,sdt)+pred_label_I][0][1] for mu in mu_list],
                 [norm_pred_I[(N,mu,nflip,sdt)+pred_label_I][2][1] for mu in mu_list],
                 c='k', ls='--', label = 'Fitness inference', lw=2)


    plt.ylabel(r'average distance $\Delta$ to future populations')
    plt.xlabel('average pairwise distance')
    #add panel label 
    plt.xscale('log')
    plt.legend(loc=1)
    plt.xlim([1.0, 100])
    plt.tight_layout()

    for ff in file_formats:
        plt.savefig(figure_folder+'Fig3_S1_pairwise_diversity_vs_distance_sdt_'+str(sdt)+'_nflip_'+str(nflip)+'_polarizer_valdt_'+str(valdt)+ff)
    #plt.close()


###################################################################
### correlation vs pairwise diversity
###################################################################


## plot gamma versus the number of predictions that are worse than random
# reload the data without averaging over the different realizations.
pred, norm_pred, run_stats,corrcoeffs = AU.load_prediction_data(prefix, N_list, mu_list, [nflip],
                                                     sdt_list, return_mean=False, polarizer=True)

for sdt in [100]:
        ### PLOT FITNESS CORRELATION VS DSCALE ###
    plt.figure(figsize= (10,6))
    ax = plt.subplot(111)
    #plt.title(r'$\omega='+str(omega)+',\;dt='+str(sdt)+'$')    
    for mi,mu in enumerate(mu_list):
        for ni, N in enumerate(N_list):
            label = (N,mu,nflip,sdt)+(valdt,)
            if ni==0:
                label_str = r'$u ='+str(mu*L)+'$'
            else:
                label_str = None
            plt.plot(m_list, [np.mean(pred[label][:,0]<pred[label][:,memi+2])
                            for memi,m  in enumerate(m_list)], lw=2, marker='o', markersize=10,
                            ls=line_styles[ni], c=cols[mi], label = label_str)

    #plt.xscale('log')
    #add panel label 
    plt.text(0.02,0.9,'Fig.~2-S3', transform = plt.gca().transAxes, fontsize = 20)
    plt.xlim([0.01, 10.5])
    plt.xscale('log')
    plt.ylabel('worse than random (out of 100)')
    plt.xlabel(r'time rescaling $\gamma$')
    plt.legend(loc=1,numpoints=1)
    for ff in file_formats:
        plt.savefig(figure_folder+'Fig2_S3_gamma_vs_predictability_polarizer_sdt_'+str(sdt)+'_nflip_'+str(nflip)+'_w_'+str(omega)+'_valdt_'+str(valdt)+ff)

