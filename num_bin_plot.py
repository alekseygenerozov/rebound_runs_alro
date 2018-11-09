from rebound_runs import bin_analysis 
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import sys
from latex_exp import latex_exp
import argparse

import h5py
from scipy.stats import sem


def num_analytic(num, v, m=5.0e-5):
	'''
	Analytic estimate for number of binaries 

	num--number of star's in sim
	v--velocity dispersion 
	m--mass of each star (5x10^-5) by default

	The disk is has an r^-3 surface density profile and extends from r=1 to r=2. (NB the corresponds to dN/da~a^-2)
	'''
	##Normalization of r^-3 surface density corresponding to a single star.
	norm=0.32
	##Evaluate vh at 1.2 to reproduce v/vh throughout the disk--more motivation?
	r1=1.2
	rh=(m/3.)**(1./3.)*r1
	vh=rh*(r1)**-1.5

	##Numerical pre-factor comes from doing integral over the disk
	return (7./8.)*(2.*np.pi/3.)/(np.pi)*num**2*norm*(4.*np.pi/3.)*rh**2.*(v/vh)**-4.

def sig_a3(r):
	return 0.32*r**-3.

def num_analytic_b(num, sig1, v_file,  m=5.0e-5):
	'''
	Applying analytic binary formula one section of the disk at a time...

	sig--function that return surface density through the disk

	'''
	tot=np.zeros(len(v_file['bin_0']))
	num_bins=len(v_file.keys())
	print num_bins
	delta_a=1./num_bins
	a_ords=np.arange(1, 2.01, delta_a)
	a_ords=np.array([0.5*(a_ords[ii]+a_ords[ii+1]) for ii in range(num_bins)])
	print a_ords
	
	for jj in range(num_bins):
		rh=((m/3.)**(1./3.)*a_ords[jj])
		vh=rh/a_ords[jj]**1.5
		# v=(v_file['bin_{0}'.format(jj)][...]['col2']**2.+v_file['bin_{0}'.format(jj)][...]['col3']**2.)**0.5
		v=v_file['bin_{0}'.format(jj)][...]['col2']

		##Number of binaries w/in each annulus
		print (2.*num*np.pi*a_ords[jj]*delta_a*sig1(a_ords[jj]))
		tot+=4.*np.pi/3.*num*sig1(a_ords[jj])*rh**2.*(v/vh)**-4*(2.*num*np.pi*a_ords[jj]*delta_a*sig1(a_ords[jj]))
	print tot
	return tot

		

parser=argparse.ArgumentParser(description='Plot number of binaries after a rebound run')
parser.add_argument('-b', '--base', help='Location of sim data')
parser.add_argument('-m', '--mass', type=float, help='Mass of each star (only used for analytic comparison)')
parser.add_argument('-y', '--ymax', type=float, default=20., help='Maximum y for plot')
parser.add_argument('-t', '--tmax', type=float, default=20., help='Maximum time for plot')
parser.add_argument('-c1', '--col1', default='black', help='Color for simulation results')
parser.add_argument('-c2', '--col2', default='red', help='Color for analytic prediction')
parser.add_argument('-na', '--nanalyt', dest='analyt',  action='store_false', help='Flag indicating whether to plot analytic solution')
parser.add_argument('-e', '--ext', default='png', help='extension for image file')


args=parser.parse_args()

base=args.base
mass=args.mass
tmax=args.tmax
ymax=args.ymax
col1=args.col1
col2=args.col2
ext=args.ext
print args.analyt
analyt=args.analyt


fig,ax=plt.subplots(figsize=(10,9))
ax.set_xlabel('Time [Orbits]')
ax.set_ylabel('Number of binaries')
ax.set_xlim(0, tmax)
ax.set_ylim(0, ymax)

t_std=np.arange(1.0e-14,(1.01)*tmax*2.*np.pi, 0.2*np.pi)
num_bins=[]
num_bins_heavy=[]
num_bins_analytic=[]
names=np.genfromtxt(base+'/names', dtype=str)
##Iterating over all runs
for ii,name in enumerate(names):
	bins=bin_analysis.BinAnalysis(base+name)
	ts=bins.ts
	if len(ts)<10.01*tmax:
		continue
	nums=bins.num_bins()

	vs=np.genfromtxt(base+name.replace('.bin', '_sigs'))
	ms=np.genfromtxt(base+name.replace('.bin', '_masses'))
	f=h5py.File(name.replace('.bin','_sigs.h5'))
	#nums_analytic = num_analytic(len(ms)-1., vs[:,2], mass)
	nums_analytic=num_analytic_b(len(ms)-1, sig_a3, f, mass)

	##Ensure number of binaries evaluate for the same grid of times
	nums=interp1d(ts, nums)(t_std)
	nums_analytic=interp1d(ts, nums_analytic[:len(ts)])(t_std)
	##Append data to list
	num_bins.append(nums)
	num_bins_analytic.append(nums_analytic)

print len(num_bins)

##Median number of binaries and standard deviation as a function of time 
nums_med=np.median(num_bins, axis=0)
err=sem(num_bins, axis=0)
ax.fill_between(t_std/(2.*np.pi), nums_med-err, nums_med+err,\
			 color=col1, alpha=0.3)
ax.plot(t_std/(2.*np.pi), nums_med, color=col1, label='Simulation')
ax.annotate('m='+'{0}'.format(latex_exp.latex_exp(mass)), (0.99*tmax,0.75*ymax), horizontalalignment='right')

##Analytic prediction
nums_med_analytic=np.median(num_bins_analytic, axis=0)
err=sem(num_bins_analytic, axis=0)
if analyt:
	ax.fill_between(t_std/(2.*np.pi), nums_med_analytic-err, nums_med_analytic+err,\
				 color=col2, alpha=0.3)
	ax.plot(t_std/(2.*np.pi), nums_med_analytic, color=col2, label='Slichting+Sari')

ax.legend()
fig.savefig(base+'/num_bins.'+ext, transparent=False)