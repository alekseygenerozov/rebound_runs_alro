from rebound_runs import bin_analysis 
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import sys
from latex_exp import latex_exp
import argparse




def num_analytic(num, v, m=5.0e-5):
	'''
	Analytic estimate for number of binaries 

	num--number of star's in sim
	v--velocity dispersion 
	m--mass of each star (5x10^-5) by default

	The disk is has an r^-3 surface density profile and extends from r=1 to r=2. (NB the corresponds to dN/da~a^-2)

	A lot of this is hard-coded for particular parameters; would be nice to generalize
	'''
	##Normalization of r^-3 surface density corresponding to a single star.
	norm=0.32
	##Measure vh at r=1.2--Not very well motivated but this gives a decent approximation v/vh throughout the disk 
	##the disk...
	r1=1.2
	rh=(m/3.)**(1./3.)*r1
	vh=rh*(r1)**-1.5

	##Numerical pre-factor comes from doing integral over the disk
	return (7./8.)*(2.*np.pi/3.)/(np.pi)*num**2*norm*(4.*np.pi/3.)*rh**2.*(v/vh)**-4.

parser=argparse.ArgumentParser(description='Plot number of binaries after a rebound run')
parser.add_argument('-b', '--base', help='Location of sim data')
parser.add_argument('-m', '--mass', type=float, help='Mass of each star (only used for analytic comparison)')
parser.add_argument('-y', '--ymax', type=float, default=20., help='Maximum y for plot')
parser.add_argument('-t', '--tmax', type=float, default=20., help='Maximum time for plot')
parser.add_argument('-c1', '--col1', default='black', help='Color for simulation results')
parser.add_argument('-c2', '--col2', default='red', help='Color for analytic prediction')
parser.add_argument('-c3', '--col3', default='green', help='Color for analytic prediction')
parser.add_argument('-mh', '--mheavy', default=1.0e-4, help='Mass separating heavies from the light')
parser.add_argument('-na', '--nanalyt', dest='analyt',  action='store_false', help='Flag indicating whether to plot analytic solution')
parser.add_argument('-e', '--ext', default='png', help='extension for image file')

args=parser.parse_args()

base=args.base
mass=args.mass
tmax=args.tmax
ymax=args.ymax
col1=args.col1
col2=args.col2
col3=args.col3
ext=args.ext
analyt=args.analyt
mheavy=args.mheavy

fig,ax=plt.subplots(figsize=(10,9))
ax.set_xlabel('Time [Orbits]')
ax.set_ylabel('Number of binaries')
ax.set_xlim(0, tmax)
ax.set_ylim(0, ymax)

t_std=np.arange(1.0e-14,(1.01)*tmax*2.*np.pi, 0.2*np.pi)
num_bins=[]
num_bins_analytic=[]
names=np.genfromtxt(base+'/names', dtype=str)
##Iterating over all runs
for ii,name in enumerate(names):
	idx=np.where(np.genfromtxt(base+name.replace('.bin', '_masses'))>mheavy)[0][0]+1
	bins=bin_analysis.BinAnalysis(base+name)
	bins_tab=bins.bins
	bins_tab_heavy=bins_tab[(bins_tab[:,1]>=idx)|(bins_tab[:,2]>=idx)]

	times_arr=bins_tab_heavy[:,0]
	nums=[len(times_arr[np.isclose(times_arr,tt, atol=0., rtol=1.0e-12)]) for tt in bins.ts]
	vs=np.genfromtxt(base+name.replace('.bin', '_sigs_heavy'))
	ms=np.genfromtxt(base+name.replace('.bin', '_masses'))

	nums_analytic = num_analytic(len(ms[ms>mheavy]), vs[:,2], mass)
	ts=bins.ts
	if len(ts)<10.01*tmax:
		continue
	##Ensure number of binaries evaluate for the same grid of times
	nums=interp1d(ts, nums)(t_std)
	nums_analytic=interp1d(ts, nums_analytic[:len(ts)])(t_std)
	##Append data to list
	num_bins.append(nums)
	num_bins_analytic.append(nums_analytic)


#Analytic prediction
nums_med_analytic=np.median(num_bins_analytic, axis=0)
err=np.std(num_bins_analytic, axis=0)
if analyt:
	ax.fill_between(t_std/(2.*np.pi), nums_med_analytic-err, nums_med_analytic+err,\
				 color=col2, alpha=0.3)
	ax.plot(t_std/(2.*np.pi), nums_med_analytic, color=col2, label='Slichting+Sari')


##Median number of binaries and standard deviation as a function of time 
nums_med=np.median(num_bins, axis=0)
err=np.std(num_bins, axis=0)
ax.fill_between(t_std/(2.*np.pi), nums_med-err, nums_med+err,\
			 color=col1, alpha=0.3)
ax.plot(t_std/(2.*np.pi), nums_med, color=col1, label='Simulation')
ax.annotate('m='+'{0}'.format(latex_exp.latex_exp(mass)), (0.99*tmax,0.75*ymax), horizontalalignment='right')


ax.legend()
fig.savefig(base+'/num_bins_heavy.'+ext, transparent=True)