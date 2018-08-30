# Created by Hayden Foote 
# JILA | Madigan Group | CU Boulder | January 2018
# Version 1.1 with one changed mass star

# This program simulates an eccentric disk of stars around a supermassive black hole

import rebound
import numpy as np
import matplotlib.pyplot as plt
import random as rand

print rebound.__version__

#---------------------------------------MAIN---------------------------------------#
bin_file=open('bins', 'ab')
from itertools import combinations
def bin_find(sim):
    # Find two closest particles
    ps = sim.particles
    m0 = ps[0].m
    for i1, i2 in combinations(range(1, sim.N),2): # get all pairs of indices/start at 1 to exclude SMBH
        dp = ps[i1] - ps[i2]   # Calculates the coponentwise difference between particles
        ##Calculate com velocity of two particles...
        ##Masses
        m1=ps[i1].m
        m2=ps[i2].m
        com = (m1*ps[i1]+m2*ps[i2])/(m1+m2)
        ##Particle velocities in com frame
        p1_com = ps[i1] - com
        p2_com = ps[i2] - com
        v12 = (p1_com.vx**2.+p1_com.vy**2.+p1_com.vz**2.)
        v22 = (p2_com.vx**2.+p2_com.vy**2.+p2_com.vz**2.)

        ##Kinetic and potential energies
        ke=0.5*m1*v12+0.5*m2*v22
        d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
        pe=sim.G*(m1*m2)/d2**0.5
        ##Hill sphere condition.
        inside_hill=((d2**0.5<(m1/m0)**(1./3.)*ps[i1].a) or (d2**0.5<(m2/m0)**(1./3.)*ps[i2].a))

        ##If the kinetic energy is less than the potential energy 
        if ((ke<pe) and (inside_hill)):
        	np.savetxt(bin_file, [[sim.t, i1, i2,  com.x, com.y, com.z,  com.vx, com.vy, com.vz]])


#Define variables
N = 100 # number of stars
pRun = 500 # number of orbital periods to run the simulation
pOut = 1 # number of orbital periods between output files

# Make a ranodmized mean anomaly array
M = np.zeros(N + 1)
for j in range(0, N + 1):
	M[j] = rand.uniform(0, 2 * np.pi)

# Make a randomized inclination array
i = np.zeros(N)
for k in range(N):
	i[k] = rand.uniform(0, 5 * np.pi / 180.0)

# Make a semimajor axis distribution
a = np.linspace(1, 2, N)

# Make the sim
sim = rebound.Simulation()
sim.G = 1.
sim.add(m = 1) # BH
# sim.add(m = 5e-3, a = 1.3, e = 0.7, inc = 0, Omega = 0, omega = 0, M = M[-1])  # One massive object, want M+Nm = .01

for l in range(0,N): # Adds stars
	sim.add(m = 5e-5, a = a[l], e = 0.7, inc = i[l], Omega = 0, omega = 0, M = M[l]) 

##Set up simulation archive...
sim.automateSimulationArchive("archive.bin",interval=2.*np.pi*pOut,deletefile=True)
# This loop runs the simulation and prints output files every pOut orbital periods
time = np.linspace(0, pRun, (pRun / pOut) + 1)

# print t 
for i in time:
	sim.move_to_com()
	sim.integrate(2 * np.pi * i)
	bin_find(sim)
	print "Reached", int(i), "orbital periods."
	
print "Done."
# Sanity check: plots end orbits
# fig = rebound.OrbitPlot(sim)
# fig.savefig("ENDtest.png")
# fig.show()
