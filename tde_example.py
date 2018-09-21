# Created by Hayden Foote/ Adapted by Aleksey Generozov
# JILA | Madigan Group | CU Boulder | January 2018

# This program simulates an eccentric disk of stars around a supermassive black hole

#import sys
#sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from rebound_runs.bin_analysis import bin_find_sim

print rebound.__version__
import sys

# Density function for semimajor axes (Hayden's implementation)
def density(min1, max1):
	xmin = 1 / max1
	xmax = 1 / min1
	x = np.linspace(xmin, xmax, num=10000)
	f = 1 / x
	rand = np.random.choice(f)
	return rand

def get_tde(sim, reb_coll):
	orbits = sim[0].calculate_orbits(primary=sim[0].particles[0])
	p1,p2 = reb_coll.p1, reb_coll.p2
	idx, idx0 = max(p1, p2), min(p1, p2)
	if idx0==0:
		##idx decremented by 1 because there is no orbit 0
		print sim[0].t, orbits[idx-1].e, idx, 'TDE!'

	return 0


#Define variables
N = 100 # number of stars
pRun = 500 # number of orbital periods to run the simulation
pOut = 1 # number of orbital periods between output files

# Make a ranodmized mean anomaly array
M = np.zeros(N + 1)
for j in range(0, N + 1):
	M[j] = rand.uniform(0, 2 * np.pi)


# Make the sim
sim = rebound.Simulation()
sim.G = 1.
sim.add(m = 1, r=1.0e-4) # BH

for l in range(0,N): # Adds stars
	sim.add(m = 5.0e-5, a = density(1.,2.), e = 0.7, inc=np.random.uniform(0, 5 * np.pi / 180.0), Omega = 0, omega = 0, M = M[l]) 

sim.collision='line'
sim.collision_resolve=get_tde

##Set up simulation archive...
sim.automateSimulationArchive("archive",interval=0.2*np.pi*pOut,deletefile=True)
# This loop runs the simulation and prints output files every pOut orbital periods

sim.move_to_com()
sim.integrate(pRun*2*np.pi)
	
print "Done."

