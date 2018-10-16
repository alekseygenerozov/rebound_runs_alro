import sys
import rebound
import numpy as np

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)
ts = [sim.t for sim in sa]
ts = np.array(ts)

np.savetxt(name.replace('.bin', '_times'), ts)