import sys
import rebound
import numpy as np

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)
ms = np.array([pp.m for pp in sa[0].particles[1:]])
ms = np.array(ms)

np.savetxt(name.replace('.bin', '_masses'), ms)