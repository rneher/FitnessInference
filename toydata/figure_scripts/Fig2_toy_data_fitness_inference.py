import glob,pickle,sys
import numpy as np
sys.path.append('../../flu/src')
import test_flu_prediction as test_flu
import matplotlib.pyplot as plt
import analysis_utils_toy_data as AU
from scipy import stats

file_formats = ['.svg', '.pdf']

# set matplotlib plotting parameters
plt.rcParams.update(test_flu.mpl_params)

line_styles = ['-', '--', '.-']
cols = ['b', 'r', 'g', 'c', 'm', 'k', 'y']
figure_folder = '../figures/'
data_dir = '../data_new'
prefix= '/20140820_'
sigma=0.03
#prefix= '/20140418_'

mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6]
N_list=[20000]
nflip_list = [0.08]
sdt_list = [1,100]
pred, norm_pred, run_stats = AU.load_prediction_data(prefix, N_list, mu_list, nflip_list,
                                                     sdt_list, return_mean=False)

mean_fitness_true_fitness_spearman_i = -4

xpos = -0.23  # panel label position
ypos = 0.95

# simulation parameters
parser = test_flu.make_toy_parser()
params = parser.parse_args()
params.pred = 'toy'
params.test = 'toy'
params.flutype='toy'
params.boost = 0.0
params.subsample = 1.0
params.gen = 12000
params.year = params.gen
params.valdt = 2
params.dt = 100
params.sample_size =200
L=2000
N=20000
nflip = 0.08
mu_sim = 32e-6
params.diffusion = 0.2
params.omega = 0.3

for gamma in [1.0]: #, 2.0, 3.0,5.0]:
    params.gamma=gamma
    for sdt in [1,100]:
        file_base = params.base_name = data_dir+'_'.join(map(str,['/N', N, 'L', L, 'nflip',nflip
                                            ,'mu',mu_sim,'sdt', sdt]))+'/'+prefix+'seqs'
        
        pred_label = params.sample_size, params.gamma, params.diffusion, params.omega, params.valdt*params.dt

        # make figure
        plt.figure(figsize= (12,5))
        ### PLOT EXAMPLE PREDICTION
        print "run example prediction"
        prediction = test_flu.predict_params(['mean_fitness'],params)
        plt.subplot(131)
        plt.text(xpos,ypos,'B', transform = plt.gca().transAxes, fontsize = 36)
        true_fitness = np.array([n.fitness for n in prediction.terminals])
        true_fitness -= true_fitness.mean()
        estimated_fitness = np.array([n.mean_fitness for n in prediction.terminals])
        estimated_fitness -= estimated_fitness.mean()
        spcorr = stats.spearmanr(true_fitness, estimated_fitness)
        im = np.zeros((prediction.nstates, prediction.nstates))
        for node, true_fit in zip(prediction.terminals, true_fitness):
            xi = np.argmin(prediction.fitness_grid<(true_fit-true_fitness.max())/sigma+prediction.fitness_grid[-2])
            im[xi,:]+=node.prob
        
        #scatter fitness
        plt.scatter(true_fitness/sigma, estimated_fitness, 
                    label = r"$\rho = "+str(np.round(spcorr[0],2))+"$") 
        plt.xlim([-3.5,2.5])
        if gamma ==0.5:
            plt.ylim([-3.5,2.5])
        else:
            plt.ylim([1.2*stats.scoreatpercentile(estimated_fitness, 2),
                      1.2*stats.scoreatpercentile(estimated_fitness, 98)])
        plt.xlabel('true fitness $[\sigma]$')
        plt.ylabel('inferred fitness $[\sigma]$')
        plt.legend(loc=4, numpoints=1)

        ### PLOT FITNESS CORRELATION VS mutation rate###
        ax = plt.subplot(132)
        muts_to_plot = mu_list[2:-1]
        sim_label = (N,mu_sim,nflip,sdt)
        plt.text(xpos,ypos,'C', transform = plt.gca().transAxes, fontsize = 36)
        plt.boxplot([pred[(N,mu,nflip,sdt)+pred_label][:,mean_fitness_true_fitness_spearman_i] for mu in muts_to_plot])
        plt.xticks(np.arange(len(muts_to_plot))+1, map(str, np.array(muts_to_plot)*L))

        plt.ylim([0,1])
        plt.ylabel(r"Spearman's correlation $\rho$")
        plt.xlabel('mutation rate [per genome]')

        ### PLOT PREDICTION SUCCESS realtive to random
        plt.subplot(133)
        plt.text(xpos,ypos,'D', transform = plt.gca().transAxes, fontsize = 36)
        plot_index = 3  # 2 = true, 3 = externa, 4 = interna;, 5 = ladder
        plt.scatter(pred[sim_label+pred_label][:,1] / pred[sim_label+pred_label][:,0],
                    pred[sim_label+pred_label][:,plot_index]/pred[sim_label+pred_label][:,0])
        print "worse than random:", np.mean(pred[sim_label+pred_label][:,plot_index]>\
                                            pred[sim_label+pred_label][:,0])
        print "within 0.1 of optimal:", np.mean(((pred[sim_label+pred_label][:,plot_index]-pred[sim_label+pred_label][:,1])/
                                                (pred[sim_label+pred_label][:,0]-pred[sim_label+pred_label][:,1]))<0.1)

        plt.plot([0,1],[0,1], c='k')
        plt.plot([0,1],[1,1], c='k', ls='--')
        plt.xlabel(r'$\Delta(\mathrm{minimal})$')
        plt.ylabel(r'$\Delta(\mathrm{prediction})$')
        plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0])
        plt.xlim([0.2,1.0])
        plt.ylim([0.2,1.5])
        plt.tight_layout()


        for ff in file_formats:
            plt.savefig(figure_folder+'Fig2_toy_data_fitness_inference_gamma_'+str(params.gamma)
                        +'_D_'+str(params.diffusion)+'_sdt_'+str(sdt)+'_valdt_'+str(params.dt*params.valdt)+ff)
        #plt.close()

##############################################################################

