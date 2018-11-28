import sys
import rebound
import numpy as np
from astropy.table import Table

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)
ms=np.genfromtxt(name.replace('.bin', '_masses'))

# sigs_light=np.empty([num_bins, len(sa), 4])
# sigs_heavy=np.empty([num_bins, len(sa), 4])
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
	Table([a1[filt_light], vs[filt_light]], names=['a','vz']).write(name.replace('.bin', '_vs.h5'), path='/bin/{0}'.format(ii), append=True, overwrite=True)
