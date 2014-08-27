import glob
import numpy as np
import matplotlib.pyplot as plt
import analysis_utils_toy_data as AU

file_formats =['.svg', '.pdf']

# set matplotlib plotting parameters
params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
plt.rcParams.update(params)


line_styles = ['-', '--', '.-']
cols = ['b', 'r', 'g', 'c', 'm', 'k', 'y']
figure_folder = '../figures/'
data_dir = '../data_new'
prefix= '/20140820_'
#prefix= '/20140418_'

N_list = [10000] #,20000]
mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6]
nflip_list = [0.01,0.02,0.04, 0.08, 0.16]
dscale_list = [0.5, 1.0, 2.0,3.0]
#nflip_list = [0.04]
sdt_list = [1,100]  #determines whether 2 genomes are sampled every generation, or 200 every 100 gen
pred, norm_pred, run_stats = AU.load_prediction_data(prefix, N_list, mu_list, nflip_list,
                                                     sdt_list, return_mean=True)

valdt = 200
ssize = 200
D = 0.5

for dscale in [0.5, 1.0, 2.0,3.0]:
    for sdt in [1,100]:
        pred_label = ssize, dscale, D, valdt

        ### PLOT FITNESS CORRELATION VS PAIRWISE DIVERSITY ###
        plt.figure(figsize= (10,6))
        ax = plt.subplot(111)

        for ni, N in enumerate(N_list):
            for fi, nflip in enumerate(nflip_list[1:]):
                plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                         [pred[(N,mu,nflip,sdt)+pred_label][0][-2] for mu in mu_list],
                         [pred[(N,mu,nflip,sdt)+pred_label][2][-2] for mu in mu_list],
                         c=cols[fi], ls=line_styles[ni], label = '$n_A = '+str(nflip)+'$', lw=2)

        plt.ylabel("Spearman's correlation")
        plt.xlabel('average pairwise distance')
        plt.xscale('log')
        plt.legend(loc=4)
        #add panel label 
        plt.text(0.02,0.9,'Fig.~2-S2', transform = plt.gca().transAxes, fontsize = 20)
        plt.xlim([0.5, 200])
        for ff in file_formats:
            plt.savefig(figure_folder+'Fig2_S2_pairwise_diversity_vs_predictability_sdt_'+str(sdt)+'_dscale_'
                    +str(dscale)+'_D_'+str(D)+'_valdt_'+str(valdt)+ff)


        ### PLOT prediction_success VS PAIRWISE DIVERSITY ###
        plt.figure(figsize= (10,6))
        ax = plt.subplot(111)

        for ni, N in enumerate(N_list):
            for fi, nflip in enumerate(nflip_list[1:]):
                plt.errorbar([run_stats[(N,mu,nflip,sdt)][-1] for mu in mu_list],
                         [norm_pred[(N,mu,nflip,sdt)+pred_label][0][1] for mu in mu_list],
                         [norm_pred[(N,mu,nflip,sdt)+pred_label][2][1] for mu in mu_list],
                         c=cols[fi], ls=line_styles[ni], label = '$n_A = '+str(nflip)+'$')

        plt.ylabel(r'distance $\bar{d}$ to future populations')
        plt.xlabel('average pairwise distance')
        #add panel label 
        plt.text(0.02,0.9,'Fig.~2-S2', transform = plt.gca().transAxes, fontsize = 20)
        plt.xscale('log')
        plt.legend(loc=1)
        plt.xlim([0.5, 200])

        for ff in file_formats:
            plt.savefig(figure_folder+'Fig2_S2_pairwise_diversity_vs_distance_sdt_'+str(sdt)+'_dscale_'+str(dscale)+'_D_'+str(D)+'_valdt_'+str(valdt)+ff)
        plt.close()


mu = mu_list[-2]
for sdt in [1,100]:
    pred_label = ssize, dscale, D, valdt

    ### PLOT FITNESS CORRELATION VS DSCALE ###
    plt.figure(figsize= (10,6))
    ax = plt.subplot(111)
    for ni, N in enumerate(N_list):
        for fi, nflip in enumerate(nflip_list[1:]):
            plt.errorbar(dscale_list,
                     [pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][0][-2] for dscale in dscale_list],
                     [pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][2][-2] for dscale in dscale_list],
                     c=cols[fi], ls='-', label = '$n_A = '+str(nflip)+'$',  lw=2)

            plt.errorbar(dscale_list,
                     [norm_pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][0][1] for dscale in dscale_list],
                     [norm_pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][2][1] for dscale in dscale_list],
                     c=cols[fi], ls='--', lw=2)

            plt.errorbar(dscale_list,
                     [norm_pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][0][2] for dscale in dscale_list],
                     [norm_pred[(N,mu,nflip,sdt)+(ssize, dscale, D, valdt)][2][2] for dscale in dscale_list],
                     c=cols[fi], ls='-.', lw=2)

    plt.text(2.05, 0.6, "inferred -- true fitness correlation")
    plt.text(2.05, 0.42, r"top ranked external nodes $\bar{d}$")
    plt.text(2.05, 0.27, r"top ranked internal nodes $\bar{d}$")
    plt.xlabel(r'scale parameter $\gamma$')
    #plt.xscale('log')
    #add panel label 
    plt.text(0.02,0.9,'Fig.~2-S3', transform = plt.gca().transAxes, fontsize = 20)
    plt.xlim([0.0, 3.5])
    plt.ylim([0.0, 0.65])
    plt.legend(loc=3)
    for ff in file_formats:
        plt.savefig(figure_folder+'Fig2_S3_dscale_vs_predictability_sdt_'+str(sdt)+'_mu_'+str(mu)+'_D_'+str(D)+'_valdt_'+str(valdt)+ff)

