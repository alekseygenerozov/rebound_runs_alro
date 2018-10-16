import sys
import rebound
import numpy as np

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)
ms=np.genfromtxt(name.replace('.bin', '_masses'))

sigs_light=np.empty([len(sa), 3])
sigs_heavy=np.empty([len(sa), 3])
for ii,sim in enumerate(sa):
	vs=np.array([pp.vxyz for pp in (sim.particles[1:])])
	filt_light=(ms<1.0e-4)
	filt_heavy=(ms>1.0e-4)
	sigs_light[ii] = np.std(vs[filt_light], axis=0)
	sigs_heavy[ii] = np.std(vs[filt_heavy], axis=0)

np.savetxt(name.replace('.bin', '_sigs_light'), sigs_light)
np.savetxt(name.replace('.bin', '_sigs_heavy'), sigs_heavy)
