#!/ebio/ag-neher/share/programs/bin/python2.7
#
#script provdies functions used to load and process the flu prediction data
#
import numpy as np
from scipy.stats import scoreatpercentile
import glob,os

def load_prediction_data_dir(dir_name, prefix,L, return_mean = True):
    predictions = {}
    normalized_predictions = {}
    if os.path.isdir(dir_name):
        file_list = glob.glob(dir_name+prefix+'seqs_prediction_*D_?.?_dt_*.dat')
        for fname in file_list:
            ssize,dscale,D,dt = map(float,fname.split('/')[-1][:-4].split('_')[4::2])
            tmp_pred = np.loadtxt(fname)[:,1:]
            min_dis = np.repeat([tmp_pred[:,1]], 4, axis=0).T
            average_dis = np.repeat([tmp_pred[:,0]], 4, axis=0).T
            #print 'min',min_dis.mean(axis=0), min_dis.shape
            #print 'avg', average_dis.mean(axis=0), average_dis.shape
            tmp_normed = (tmp_pred[:,2:-2]-min_dis)/(average_dis-min_dis+0.01/L)
            #print 'norm',tmp_normed.mean(axis=0)
            if return_mean:
                predictions[(ssize,dscale, D,dt)] = (tmp_pred.mean(axis=0), tmp_pred.std(axis=0), tmp_pred.std(axis=0)/np.sqrt(tmp_pred.shape[0]))
                normalized_predictions[(ssize,dscale, D,dt)] = (tmp_normed.mean(axis=0),
                            tmp_normed.std(axis=0), tmp_normed.std(axis=0)/np.sqrt(tmp_pred.shape[0]))
            else:
                predictions[(ssize,dscale, D,dt)] = tmp_pred
                normalized_predictions[(ssize,dscale, D,dt)] = tmp_normed
        run_stats = np.loadtxt(dir_name+prefix+'stats.dat')[500:,1:].mean(axis=0)
        return predictions, normalized_predictions,run_stats
    else:
        print dir_name,'does not exist'
        return None,None,None
    
def load_prediction_data(prefix, N_list = [10000, 20000], mu_list = [1e-6,2e-6, 4e-6, 8e-6, 16e-6, 32e-6, 64e-6, 128e-6],
                         nflip_list = [0.01,0.02,0.04, 0.08, 0.16], sdt_list = [100], L=2000, return_mean = True):
    predictions = {}
    run_stats = {}
    normalized_predictions= {}
    for N in N_list:
        for mu in mu_list:
            for nflip in nflip_list:
                for sdt in sdt_list:
                    dirname = '_'.join(map(str,['../data_new/N', N, 'L', L, 'nflip',nflip,'mu',mu,'sdt', sdt]))
                    tmp_pred, tmp_normed, tmp_run_stats = load_prediction_data_dir(dirname, prefix,L, return_mean)
                    if tmp_pred is not None:
                        for p,val in tmp_pred.iteritems():
                            predictions[(N,mu,nflip,sdt)+p] = val
                        for p,val in tmp_normed.iteritems():
                            normalized_predictions[(N,mu,nflip,sdt)+p] = val
                        run_stats[(N,mu,nflip,sdt)] = tmp_run_stats
    return predictions, normalized_predictions, run_stats
    
