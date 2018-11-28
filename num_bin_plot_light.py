from rebound_runs import bin_analysis 
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import sys
from latex_exp import latex_exp
import argparse

from extrap import extrap
from scipy.stats import sem

from scipy.special import erf

def num_analytic_gen(num, v, mu, s, m=5.0e-5):
	'''
	Analytic estimate for number of binaries 

	num--number of star's in sim
	v--velocity dispersion 
	m--mass of each star (5x10^-5) by default

	'''
	r1=1.2
	rh0=(m/3.)**(1./3.)
	rh=rh0*r1
	vh=rh*(r1)**-1.5

	##Numerical pre-factor comes from doing integral of sigma^2* over the disk
	#integ=(np.exp((9*mu**2)/4. - 3*s)*np.sqrt(mu**(-2))*rh**2*(1 + erf((np.sqrt(mu**(-2))*(3*mu**2 - 2*s))/2.)))/(4.*np.sqrt(np.pi))
	integ=(np.exp(mu**2/4. - s))/(8.*mu*np.pi**2.5)

	return (4.*np.pi/3.)*integ*num**2*rh0**2.*(v/vh)**-4.


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
	r1=1.2
	rh0=(m/3.)**(1./3.)
	rh=rh0*r1	
	vh=rh*(r1)**-1.5

	##Numerical pre-factor comes from doing integral over the disk--Should in principlel replace rh with rh0, but
	##this makes the agreement worse not better.  
	return (7./8.)*(2.*np.pi/3.)/(np.pi)*num**2*norm*(4.*np.pi/3.)*rh**2.*(v/vh)**-4.

def num_analytic_b(num, vratio, m=5.0e-5):
	'''
	Analytic estimate for number of binaries 

	num--number of star's in sim
	v--velocity dispersion 
	m--mass of each star (5x10^-5) by default

	The disk is has an r^-3 surface density profile and extends from r=1 to r=2. (NB the corresponds to dN/da~a^-2)
	'''
	##Normalization of r^-3 surface density corresponding to a single star.
	norm=0.32
	r1=1.2
	rh=(m/3.)**(1./3.)*r1
	# vh=rh*(r1)**-1.5

	##Numerical pre-factor comes from doing integral over the disk
	return (7./8.)*(2.*np.pi/3.)/(np.pi)*num**2*norm*(4.*np.pi/3.)*rh**2.*(vratio)**-4.

parser=argparse.ArgumentParser(description='Plot number of binaries after a rebound run')
parser.add_argument('-b', '--base', help='Location of sim data')
parser.add_argument('-m', '--mass', type=float, help='Mass of each star (only used for analytic comparison)')
parser.add_argument('--ymax', type=float, default=20., help='Maximum y for plot')
parser.add_argument('--ymin', type=float, default=20., help='Maximum y for plot')
parser.add_argument('-t', '--tmax', type=float, default=20., help='Maximum time for plot')
parser.add_argument('-c1', '--col1', default='black', help='Color for simulation results')
parser.add_argument('-c2', '--col2', default='red', help='Color for analytic prediction')
parser.add_argument('-c3', '--col3', default='green', help='Color for analytic prediction')
# parser.add_argument('-mh', '--mheavy', default=1.0e-4, help='Mass separating heavies from the light')
parser.add_argument('-e', '--ext', default='png', help='extension for image file')



args=parser.parse_args()

base=args.base
mass=args.mass
tmax=args.tmax
ymax=args.ymax
ymin=args.ymin
col1=args.col1
col2=args.col2
col3=args.col3
ext=args.ext
# mheavy=args.mheavy

fig,ax=plt.subplots(figsize=(10,9))
ax.set_xlabel('Time [Orbits]')
ax.set_ylabel('Number of binaries')
ax.set_xlim(1, tmax)
ax.set_ylim(ymin, ymax)

t_std=np.arange(1.0e-14,(1.01)*tmax*2.*np.pi, 0.2*np.pi)
num_bins=[]
# num_bins_analytic=[]
names=np.genfromtxt(base+'/names', dtype=str)
vs_all=[]
nstars=[]
##Iterating over all runs
for ii,name in enumerate(names):
	bins=bin_analysis.BinAnalysis(base+name)
	bins_tab=bins.bins
	vs=np.genfromtxt(base+name.replace('.bin', '_sigs_low'))
	##Taking only the z-component of the velocity dispersion
	vs=vs[:,2]
	ts=bins.ts

	ms=np.genfromtxt(base+name.replace('.bin', '_masses'))
	mheavy=np.median(ms)

	##Select only the light-light binaries
	idx=np.where(ms<=mheavy)[0][-1]+1
	bins_tab_light=bins_tab[(bins_tab[:,1]<idx) & (bins_tab[:,2]<idx)]
	times_arr=bins_tab_light[:,0]
	nstars.append(idx)
	
	# try:
	# 	v_arr=extrap.extrap1d(interp1d(ts, vs))(times_arr)
	# except Exception as e:
	# 	print e.message
	# 	continue
	# a_thres=2.*mass/3./v_arr**2.
	# # print a_thres[-1]
	# filt_a=bins_tab_light[:,4]>=a_thres
	# bins_tab_light=bins_tab_light[filt_a]
	# times_arr=bins_tab_light[:,0]



	nums=[len(times_arr[np.isclose(times_arr,tt, atol=0., rtol=1.0e-12)]) for tt in ts]
	#vvh=np.genfromtxt(name.replace('.bin', '_vvh_ratio_low'))
	# nums_analytic = num_analytic(len(ms[ms<=mheavy]), vs, mass)
	#nums_analytic = num_analytic_b(len(ms[ms<=mheavy]), vvh, mass)
	if len(ts)<10.01*tmax:
		continue
	vs_all.append(interp1d(ts, vs)(t_std))
	##Ensure number of binaries evaluate for the same grid of times
	nums=interp1d(ts, nums)(t_std)
	# nums_analytic=interp1d(ts, nums_analytic[:len(ts)])(t_std)
	##Append data to list
	num_bins.append(nums)
	# num_bins_analytic.append(nums_analytic)

vs=np.mean(vs_all, axis=0)
print np.mean(nstars)
nums_med_analytic=num_analytic(np.mean(nstars), vs, mass)
print nums_med_analytic


# print nums_med_analytic[-1]
# err=sem(num_bins_analytic, axis=0)
# ax.fill_between(t_std/(2.*np.pi), nums_med_analytic-err, nums_med_analytic+err,\
# 			 color=col2, alpha=0.3)
ax.loglog(t_std/(2.*np.pi), nums_med_analytic, color=col2, label='Slichting+Sari')


##Median number of binaries and standard deviation as a function of time 
nums_med=np.mean(num_bins, axis=0)
err=sem(num_bins, axis=0)
ax.fill_between(t_std/(2.*np.pi), nums_med-err, nums_med+err,\
			 color=col1, alpha=0.3)
ax.loglog(t_std/(2.*np.pi), nums_med, color=col1, label='Simulation')
# ax.annotate('m='+'{0}'.format(latex_exp.latex_exp(mass)), (0.99*tmax,0.75*ymax), horizontalalignment='right')


ax.legend()
plt.show()
# fig.savefig(base+'/num_bins_light.'+ext)