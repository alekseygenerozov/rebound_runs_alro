# Created by Hayden Foote/ Adapted by Aleksey Generozov
# JILA | Madigan Group | CU Boulder | January 2018

# This program simulates an eccentric disk of stars around a supermassive black hole


import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from rebound_runs.bin_analysis import bin_find_sim

print rebound.__version__

# Density function for semimajor axes (Hayden's implementation)
def density(min1, max1):
    xmin = 1 / max1
    xmax = 1 / min1
    x = np.linspace(xmin, xmax, num=10000)
    f = 1 / x
    rand = np.random.choice(f)
    return rand

# Aleksey's implementation.
def density2(min1, max1):
    r=np.random.random(1)[0]
    return (1./min1-r*(1./min1-1./max1))**-1.

#Define variables
N = 110 # number of stars
pRun = 500 # number of orbital periods to run the simulation
pOut = 1 # number of orbital periods between output files

# Make a ranodmized mean anomaly array
M = np.zeros(N + 1)
for j in range(0, N + 1):
	M[j] = rand.uniform(0, 2 * np.pi)


# Make the sim
sim = rebound.Simulation()
sim.G = 1.
sim.add(m = 1) # BH
# sim.add(m = 5e-3, a = 1.3, e = 0.7, inc = 0, Omega = 0, omega = 0, M = M[-1])  # One massive object, want M+Nm = .01

for l in range(0,N): # Adds stars
	sim.add(m = 5e-5, a = density(1.,2.), e = 0.7, inc=np.random.uniform(0, 5 * np.pi / 180.0), Omega = 0, omega = 0, M = M[l]) 

##Delete primordial binaries
bins=bin_find_sim(sim)
bins=np.array(bins)
##Delete in reverse order (else the indices would become messed up)
to_del=np.unique((np.sort(bins[:,1])[::-1]).astype(int))
print to_del
print len(sim.particles)
for idx in to_del:
    sim.remove(idx)
print len(sim.particles)

##Set up simulation archive...
sim.automateSimulationArchive("archive.bin",interval=0.2*np.pi*pOut,deletefile=True)
# This loop runs the simulation and prints output files every pOut orbital periods

sim.move_to_com()
sim.integrate(pRun*2*np.pi)
	
print "Done."
# Sanity check: plots end orbits
# fig = rebound.OrbitPlot(sim)
# fig.savefig("ENDtest.png")
# fig.show()
