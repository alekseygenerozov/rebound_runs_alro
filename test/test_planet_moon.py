from bash_command import bash_command as bc

import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
from rebound_runs import bin_analysis
import numpy as np



a_bin=3.0e-3
inc_bin=np.pi/3.
e_bin=0.1
tol=1.0e-5

sim2 = rebound.Simulation()
sim2.G = 1.

sim2.add(m = 1.) # Star
sim2.add(m = 3e-6, a=1., inc=0.) #Planet
sim2.add(m = 3e-10, a=a_bin, e=e_bin, inc=inc_bin, primary=sim2.particles[1]) #Moon

bc.bash_command('rm test.bin')
sim2.simulationarchive_snapshot('test.bin')


def test_bin_props():
	x=bin_analysis.bin_props(sim2.particles[1], sim2.particles[2])
	assert abs(x[1]-a_bin)/a_bin<tol
	assert abs(x[2]**0.5-e_bin)/e_bin<tol
	print x[-1]*np.pi/180.
	assert abs(x[-1]*np.pi/180.-inc_bin)/inc_bin<tol

def test_bin_find():
	bins=bin_analysis.bin_find([0, 'test.bin'])
	print bins
	assert (len(bins)==1)
	assert set(bins[0,[1,2]])=={1,2}


