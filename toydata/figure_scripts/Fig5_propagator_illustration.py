import sys
sys.path.append('../../prediction_src/')
import solve_survival
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10

file_formats = ['.svg','.pdf']

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

D=0.2
eps=0.001
t_end = 6
nt = 20
fitness_grid = np.linspace(-4,4,81)
t = np.linspace(0,t_end-0.1,nt)
mysol = solve_survival.survival_gen_func(fitness_grid)
mysol.integrate_phi(D,eps,  np.linspace(0,20,101),dt=0.005)

# solve for the propagator with several starting (closer to present) times, always the same endtime
# save results in arrays for backward and forward separately
x = mysol.fitness_grid[1:-1]
Z_1 = np.zeros((nt, x.shape[0]))
Z_2 = np.zeros((nt, x.shape[0]))
Y_1 = np.zeros((nt, mysol.fitness_grid.shape[0]))
Y_2 = np.zeros((nt, mysol.fitness_grid.shape[0]))
y1 = fitness_grid[0.5*fitness_grid.shape[0]-1]  # y1, y2 are the initial fitness of the ancestor
y2 = fitness_grid[0.8*fitness_grid.shape[0]-1]  # -1 since the only 1:-1 of the fitness grid is used
for ti,t_start in enumerate(t):
    propagator = mysol.integrate_prop(D,eps, mysol.fitness_grid[4:-4], t_start, t_end,dt=0.005)
    Z_1[ti] = propagator[-1,:,0.5*fitness_grid.shape[0]]
    Z_2[ti] = propagator[-1,:,0.8*fitness_grid.shape[0]]
    Y_1[ti] = propagator[-1,0.2*fitness_grid.shape[0],:]
    Y_2[ti] = propagator[-1,0.8*fitness_grid.shape[0],:]
    Y_1[ti]/=Y_1[ti].sum()*mysol.dx
    Y_2[ti]/=Y_2[ti].sum()*mysol.dx

# make propagator figure. 3 panels.  
fig= plt.figure(figsize = (12,5))
xpos = -0.23  # panel label position
ypos = 0.95


# Panel A: distribution of offspring
ax=plt.subplot(131)
plt.text(xpos,ypos,'A', transform = plt.gca().transAxes, fontsize = 36)
#plt.title('offspring fitness')
for ti in range(len(t)-1,11,-1):
    if ti%1==0:
        plt.plot(x, Z_2[ti,:], 
                c = mpl.cm.Reds((ti-10)*255/10), lw=3, alpha=0.8)
for ti in range(len(t)-1,11,-1):
    if ti%1==0:
        plt.plot(x, Z_1[ti,:], 
                c = mpl.cm.Blues((ti-10)*255/10), lw=3, alpha=0.8)

plt.yscale('log')
plt.ylim([0.01,30])
plt.ylabel('offspring fitness distribution')
plt.xlabel('fitness $[\sigma]$')

## panel B: reproductive value
ax=plt.subplot(132)
plt.text(xpos,ypos,'B', transform = plt.gca().transAxes, fontsize = 36)
#plt.title('offspring number')
rep_val_2 = [ np.sum(Z_2[ti,:])*mysol.dx for ti in range(len(t))[::-1]]
plt.plot(t, rep_val_2, c='r')
plt.plot(t, np.exp(y2*t - t**2/2 + D*t**3/3), ls='--', c='r')
rep_val_1 = [ np.sum(Z_1[ti,:])*mysol.dx for ti in range(len(t))[::-1]]
plt.plot(t, rep_val_1, c='b')
plt.plot(t, np.exp(y1*t - t**2/2 + D*t**3/3), ls='--', c='b')
plt.yscale('log')
plt.xlim([0,5])
plt.ylim([3e-3,3e2])
plt.ylabel('reproductive value')
plt.xlabel('time $[\sigma^{-1}]$')


ax=plt.subplot(133)
plt.text(xpos,ypos,'C', transform = plt.gca().transAxes, fontsize = 36)

#plt.title('ancestor fitness')
for ti in range(len(t)-1,-1,-1):
    if ti%2==0:
        plt.plot(mysol.fitness_grid, Y_1[ti,:], 
                c = mpl.cm.Blues(ti*255/nt), lw=3, alpha=0.8)
for ti in range(len(t)-1,-1,-1):
    if ti%2==0:
        plt.plot(mysol.fitness_grid, Y_2[ti,:], 
                c = mpl.cm.Reds(ti*255/nt), lw=3, alpha=0.8)

plt.yscale('log')
plt.ylim([0.001,10])
plt.ylabel('ancestral fitness distribution')
plt.xlabel('fitness $[\sigma]$')
plt.tight_layout()
#plt.show()
for ff in file_formats:
    plt.savefig('../figures/Fig5_propagator_illustration'+ff)



