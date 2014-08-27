import glob,pickle,sys
import numpy as np
sys.path.append('../../flu/src')
import test_flu_prediction as flu
import matplotlib.pyplot as plt
import analysis_utils_toy_data as AU
from scipy import stats

file_formats = ['.svg', '.pdf']

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
sigma=0.03
#prefix= '/20140418_'

N_list = [10000]
mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6]
#nflip_list = [0.01,0.02,0.04, 0.08, 0.16]
nflip_list = [0.04]
sdt_list = [1,100]
pred, norm_pred, run_stats = AU.load_prediction_data(prefix, N_list, mu_list, nflip_list,
                                                     sdt_list, return_mean=False)


xpos = -0.23  # panel label position
ypos = 0.95

# simulation parameters
N=10000
L=2000
nflip = 0.04
dt=100
year = 9000
mu_sim = 6.4e-5
valdt = 200
ssize = 200
D = 0.5

for dscale in [0.5, 1.0, 2.0, 3.0]:
    for sdt in [1,100]:
        
        pred_label = ssize, dscale, D,valdt

        # make figure
        plt.figure(figsize= (12,5))
        ### PLOT EXAMPLE PREDICTION
        print "run example prediction"
        file_base = data_dir+'_'.join(map(str,['/N', N, 'L', L, 'nflip',nflip,'mu',mu_sim,'sdt', sdt]))+\
            '/'+prefix+'seqs'
        with open(file_base+'_annotation.pickle', 'r') as infile:
            annotation = pickle.load(infile)

        with open(file_base+'_outgroup.fa', 'r') as infile:
            from Bio import SeqIO
            outgroup = SeqIO.read(infile,'fasta')

        aln_fname = file_base+'_AT.fasta.gz'
        cds = None
        prediction_set={'start':year-0.5*dt, 'stop':year+0.5*dt, 'regions': ['toy'], 'sample_size':ssize}

        prediction = flu.predict(aln_fname, outgroup, annotation,\
                                ['mean_fitness'],
                                prediction_set, cds, time_bins = None, subsample_factor = 1.0, boost = 0.0,
                                eps_branch_length = 1e-5, collapse = False, dscale = dscale, D=D)
        plt.subplot(131)
        plt.text(xpos,ypos,'B', transform = plt.gca().transAxes, fontsize = 36)
        true_fitness = np.array([n.fitness for n in prediction.terminals])
        true_fitness -= true_fitness.mean()
        estimated_fitness = np.array([n.mean_fitness for n in prediction.terminals])
        estimated_fitness -= estimated_fitness.mean()
        spcorr = stats.spearmanr(true_fitness, estimated_fitness)

        #scatter fitness
        plt.scatter(true_fitness/sigma, estimated_fitness, 
                    label = r"$\rho = "+str(np.round(spcorr[0],2))+"$") 
        plt.xlim([-3.5,2.5])
        if dscale ==0.5:
            plt.ylim([-3.5,2.5])
        else:
            plt.ylim([1.2*stats.scoreatpercentile(estimated_fitness, 5),
                      1.2*stats.scoreatpercentile(estimated_fitness, 95)])
        plt.xlabel('true fitness $[\sigma]$')
        plt.ylabel('inferred fitness $[\sigma]$')
        plt.legend(loc=4, numpoints=1)

        ### PLOT FITNESS CORRELATION VS mutation rate###
        ax = plt.subplot(132)
        plt.text(xpos,ypos,'C', transform = plt.gca().transAxes, fontsize = 36)
        plt.boxplot([pred[(N,mu,nflip,sdt)+pred_label][:,-2] for mu in mu_list[3:]])
        plt.xticks(np.arange(len(mu_list[3:]))+1, map(str, np.array(mu_list[3:])*L))

        plt.ylim([0,1])
        plt.ylabel(r"Spearman's correlation $\rho$")
        plt.xlabel('mutation rate [per genome]')

        ### PLOT PREDICTION SUCCESS realtive to random
        plt.subplot(133)
        plt.text(xpos,ypos,'D', transform = plt.gca().transAxes, fontsize = 36)
        plot_index = 3  # 2 = true, 3 = externa, 4 = interna;, 5 = ladder
        sim_label = (N,mu_sim,nflip,sdt)
        plt.scatter(pred[sim_label+pred_label][:,1] / pred[sim_label+pred_label][:,0],
                    pred[sim_label+pred_label][:,plot_index]/ pred[sim_label+pred_label][:,0])
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
            plt.savefig(figure_folder+'Fig2_toy_data_fitness_inference_dscale_'+str(dscale)+'_D_'+str(D)+'_sdt_'+str(sdt)+'_valdt_'+str(valdt)+ff)
        plt.close()

##############################################################################

