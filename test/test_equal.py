from bash_command import bash_command as bc

import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
from rebound_runs import bin_analysis
import numpy as np

sim2 = rebound.Simulation()
sim2.G = 1.
a_bin=1.
e_bin=0.5
tol=1.0e-5

sim2.add(m = 1.) # Star
sim2.add(m = 1, a=1, e=0.5) #Star

def test_equal():
	x=bin_analysis.bin_props(sim2.particles[0], sim2.particles[1])
	print x
	assert abs(x[1]-a_bin)/a_bin<tol
	assert abs(x[2]**0.5-e_bin)/e_bin<tol