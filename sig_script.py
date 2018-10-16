import sys
import rebound
import numpy as np

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)
sigs = [np.std([np.array(pp.vxyz) for pp in sim.particles[1:]],axis=0) for sim in sa]
sigs = np.array(sigs)

np.savetxt(name.replace('.bin', '_sigs'), sigs)