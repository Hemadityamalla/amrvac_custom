#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
from subprocess import check_output
from hemaditya_euler import *

p = argparse.ArgumentParser()
p.add_argument('-scheme', type=str,
               choices=['tvdlf', 'tvdmu', 'hll', 'hllc'],
               default='tvdlf', help='Which scheme to use')
p.add_argument('-limiter', type=str,
               choices=['minmod', 'woodward', 'mcbeta', 'superbee', 'albada',
                        'koren', 'vanleer', 'cada', 'cada3', 'mp5'],
               default='minmod', help='Which limiter to use')
p.add_argument('-integrator', type=str,
               choices=['twostep', 'twostep_trapezoidal', 'threestep',
                        'fourstep', 'ssprk54', 'ssprk43', 'rk4', 'jameson'],
               default='threestep', help='Which integrator to use')
p.add_argument('-tfinal', type=float, default=0.15,
               help='final time and the time at which the plot is generated')
p.add_argument('-nx', type=int, default=128,
               help='Number of cells')
               
p.add_argument('-prob', type=int, default=1, help='Which initial conditions')
args = p.parse_args()


with open('this_scheme.par', 'w') as f:
    f.write('! Temporary par file\n')
    f.write('&savelist\n')
    f.write(' dtsave_log = {}\n'.format(args.tfinal/10.0))
    f.write(' dtsave_dat = {}\n'.format(args.tfinal/10.0))
    f.write('/\n')
    f.write('&stoplist\n')
    f.write(' time_max = {}\n'.format(args.tfinal))
    f.write('/\n')
    f.write('&methodlist\n')
    f.write(' time_integrator = "{}"\n'.format(args.integrator))
    f.write(' flux_scheme = "{}"\n'.format(args.scheme))
    f.write(' limiter = "{}"\n'.format(args.limiter))
    f.write('/\n')
    f.write('&meshlist\n')
    f.write(' domain_nx1 = {}\n'.format(args.nx))
    f.write(' iprob = {}\n'.format(args.prob))
    f.write('/\n')

print('Running AMRVAC code..\n')
res = check_output(["./amrvac", "-i", "sod.par", "this_scheme.par"])

print('Running Hemadityas code..\n')
#Initial conditions
initconds = {1:np.array([[1.0,0.0,1.0],[0.125,0.0,0.1]]),\
	     2:np.array([[0.445,0.698,3.527],[0.5,0.0,0.571]]),\
	     3:np.array([[1.0,-0.5,1.0],[1.0,0.5,1.0]]),\
	     4:np.array([[1.0,0.0,1000],[1.0,0.0,0.01]]),\
	     5:np.array([[1.0,0.0,0.01],[1.0,0.0,100.0]]),\
	     6:np.array([[5.9924,19.5975,460.894],[5.99242,-6.19633,46.095]])
	    }

ics = initconds.get(args.prob)
gamma = 1.4
def returnInternalEnergy(data):
	return data[:,2]/(data[:,0]*(gamma - 1.0))
	

hema = euler(args.nx, 0.3, ics, args.tfinal, gamma)
data = np.genfromtxt('output/sod0010.blk', skip_header=3)

xvals = np.copy(data[:,0])
data = np.c_[data, returnInternalEnergy(data[:,1:])]
hema[1] = np.c_[hema[1], returnInternalEnergy(hema[1])]

fig, ax = plt.subplots(2,2, figsize=(12,12))
names = ['Density','Velocity','Pressure','Sp. internal energy']; symbol='.'
#marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=15, markerfacecoloralt='tab:red')
#ax[0,0].plot(xvals, data[:,1],'r'+symbol, hema[0], hema[1][:,0]), ax[0,0].set_xlabel('x'), ax[0,0].set_ylabel(names[0])
#ax[0,1].plot(xvals, data[:,2],'b'+symbol, hema[0], hema[1][:,1]), ax[0,1].set_xlabel('x'), ax[0,1].set_ylabel(names[1])
#ax[1,0].plot(xvals, data[:,3],'y'+symbol, hema[0], hema[1][:,2]), ax[1,0].set_xlabel('x'), ax[1,0].set_ylabel(names[2])
#ax[1,1].plot(xvals, returnInternalEnergy(data[:,1:]),'y'+symbol, hema[0], returnInternalEnergy(hema[1])), ax[1,1].set_xlabel('x'), ax[1,1].set_ylabel(names[3])

for i in range(2):
	for j in range(2):
		linIdx = 2*i + j
		#print(i,j)
		ax[i,j].plot(xvals, data[:,linIdx+1],'r'+symbol, hema[0], hema[1][:,linIdx]), ax[i,j].set_xlabel('x'), ax[i,j].set_ylabel(names[linIdx])
plt.show()	
'''
plt.plot(data[:, 0], data[:, 1],'o', hema[0], hema[1][:,0])
plt.xlabel('x'); plt.ylabel('Density')
plt.subplot(222)
plt.plot(data[:, 0], data[:, 2],'o', hema[0], hema[1][:,1])
plt.subplot(223)
plt.plot(data[:, 0], data[:, 3],'o', hema[0], hema[1][:,2])
plt.subplot(224)
plt.plot(data[:, 0], returnInternalEnergy(data[:,1:]),'o', hema[0], returnInternalEnergy(hema[1]))
plt.show()
'''











