from rebound_runs import bin_analysis
import rebound
import sys
import numpy as np

from astropy.table import Table

name=sys.argv[1]
bins=bin_analysis.BinAnalysis(name)
##Set up simulation archive
#--------------------------------------------------------------------------------------------------_#
sa=rebound.SimulationArchive(name)
# ms = np.array([pp.m for pp in sa[0].particles[1:]])
ms=np.genfromtxt(name.replace('.bin', '_masses'))
#--------------------------------------------------------------------------------------------------_#
##Velocity dispersions to file
sigs=np.empty([len(sa), 3])
sigs_low=np.empty([len(sa), 3])
sigs_high=np.empty([len(sa), 3])
for ii, sim in enumerate(sa):
	vs = np.array([np.array(pp.vxyz) for pp in sim.particles[1:]])
	sigs[ii] = np.std(vs,axis=0) 
	np.savetxt(name.replace('.bin', '_sigs'), sigs)

	sigs_high[ii] = np.std(vs[ms>np.median(ms)], axis=0)
	sigs_low[ii] = np.std(vs[ms<=np.median(ms)], axis=0)
	
np.savetxt(name.replace('.bin', '_sigs'), sigs)
np.savetxt(name.replace('.bin', '_sigs_high'), sigs_high)
np.savetxt(name.replace('.bin', '_sigs_low'), sigs_low)
##Save orbital elements to hdf5 file.
#--------------------------------------------------------------------------------------------------_#
elem_names=['a', 'e', 'inc', 'omega']
#elem_names=['a', 'e']
interval=1
for ii in range(0, len(sa), interval):
	sim=sa[ii]
	nn=len(sim.particles)

	orbits=sim.calculate_orbits(primary=sim.particles[0])
	elems=np.array([getattr(orbits[jj],en) for jj in range(nn-1) for en in elem_names]).reshape([nn-1, len(elem_names)])
	tab=Table(elems, names=elem_names)	
	
	tab.write(name.replace('.bin','_elems.hdf5'), '/{0}'.format(ii), format='hdf5', append=True, overwrite=True)
