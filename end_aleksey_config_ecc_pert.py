#!/usr/bin/env python

import ConfigParser
import argparse
import uuid
# import sys

import numpy as np
from collections import OrderedDict
# sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import random as rand
from rebound_runs.bin_analysis import bin_find_sim

# Density function for semimajor axes (Hayden's implementation)
def density(min1, max1):
	xmin = 1 / max1
	xmax = 1 / min1
	x = np.linspace(xmin, xmax, num=10000)
	f = 1 / x
	rand = np.random.choice(f)
	return rand

# Aleksey's implementation.
# def density2(min1, max1):
#     r=np.random.random(1)[0]
#     return (1./min1-r*(1./min1-1./max1))**-1.

def heartbeat(sim):
# sim is a pointer to the simulation object,
# thus use contents to access object data.
# See ctypes documentation for details.
	print(sim.contents.dt)

def main():
	parser=argparse.ArgumentParser(
		description='Set up a rebound run')
	parser.add_argument('config', nargs=1,
		help='File containing config information for our run.')

	args=parser.parse_args()
	config_file=args.config[0]
	tag=str(uuid.uuid4())

	config=ConfigParser.SafeConfigParser(defaults={'name': 'archive'.format(tag), 'N':'100', 'e':'0.7',\
		'a_min':'1.', 'a_max':'2.', 'm':'5e-5', 'pRun':'500', 'pOut':'0.2'}, dict_type=OrderedDict)
	# config.optionxform=str
	config.read(config_file)

	name=config.get('params', 'name')
	name=name+"_"+tag+".bin"
	N=int(config.get('params', 'N'))
	e=config.getfloat('params', 'e')
	m=config.getfloat('params', 'm')
	a_min=config.getfloat('params', 'a_min')
	a_max=config.getfloat('params', 'a_max')
	pRun=config.getfloat('params', 'pRun')
	pOut=config.getfloat('params', 'pOut')
	print name, N, e, m, a_max, a_min, pOut

	M = np.zeros(N + 1)
	for j in range(0, N + 1):
		M[j] = rand.uniform(0, 2 * np.pi)

	sim = rebound.Simulation()
	sim.G = 1.	
	sim.add(m = 1) # BH

	for l in range(0,N): # Adds stars
		sim.add(m = m, a = density(a_min, a_max), e = rand.uniform(e-0.05, e+0.05), inc=np.random.uniform(0, 5 * np.pi / 180.0), Omega = 0, omega = 0, M = M[l])

	##Integrate forward a small amount time to initialize accelerations.
	sim.integrate(1.0e-10)
	bins=bin_find_sim(sim)
	bins=np.array(bins)
	##Delete in reverse order (else the indices would become messed up)
	if len(bins>0):
		to_del=(np.sort(np.unique(bins[:,1]))[::-1]).astype(int)
		print len(to_del)
		print len(sim.particles)
		for idx in to_del:
			sim.remove(idx)
		print len(sim.particles) 

	sim.automateSimulationArchive(name,interval=np.pi*pOut,deletefile=True)
	sim.heartbeat=heartbeat
	sim.move_to_com()
	print rebound.__version__
	sim.integrate(pRun*2*np.pi)




if __name__ == '__main__':
	main()








