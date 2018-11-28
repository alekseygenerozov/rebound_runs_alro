import sys
import rebound
import numpy as np
from astropy.table import Table

name=sys.argv[1]
num_bins=sys.argv[2]
num_bins=int(num_bins)
print name
sa=rebound.SimulationArchive(name)
ms=np.genfromtxt(name.replace('.bin', '_masses'))

sigs_light=np.empty([num_bins, len(sa), 4])
sigs_heavy=np.empty([num_bins, len(sa), 4])
for ii,sim in enumerate(sa):
	print ii
	##Orbital smas and eccentricities
	orbits = sim.calculate_orbits(primary=sim.particles[0])
	odat=np.array([[oo.a, oo.e] for oo in orbits])
	a1=odat[:,0]
	eccs=odat[:,1]
	##Look at just the light stars
	filt_light=(ms<=np.median(ms))
	# filt_heavy=(ms>1.0e-4)

	vs=np.array([pp.vz for pp in (sim.particles[1:])])
	for jj in range(num_bins):
		##Bin particles by sma and compute the velocity in each bin -- May be good to try binning by radius as well 
		a_in=1.+(1./num_bins)*jj
		a_out=1.+(1./num_bins)*(jj+1.)

		filt=((a1>=a_in) & (a1<a_out))
		# if len(eccs[filt & filt_light])==0:
		# 	tmp=np.inf
		# 	sigs_light[jj,ii]=np.array([np.inf, np.inf, np.inf, np.inf])
		# 	continue
		# tmp = np.mean(eccs[filt & filt_light]**2.)**0.5/2.**0.5*(1./(np.mean([a_in, a_out])))**0.5
		# sigs_light[jj,ii] = np.hstack([np.std(vs[filt & filt_light], axis=0), tmp])
		Table([a1[filt_light], vs[filt_light]], names=['a','vz']).write(name.replace('.bin', '_vs.h5'), path='/bin/{0}'.format(ii), append=True, overwrite=True)
		#np.savetxt(name.replace('.bin', '_vs_{0}_{1}'.format(ii, jj)), vs[filt])