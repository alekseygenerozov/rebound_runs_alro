from bash_command import bash_command as bc
import rebound
import numpy as np
from scipy.stats import kstest


@np.vectorize
def cum_dist(x):
    if x<1:
        return 0.
    elif x>=3:
        return 1.
    else:
        return (1.-1./x)/(1.-1./2.)

def test_init():
	bc.bash_command('rm archive*')
	bc.bash_command('python ../end_aleksey_config_b.py')
	name=bc.bash_command("echo archive*bin").replace('\n', '')
	sa=rebound.SimulationArchive(name)
	sim=sa[0]
	##100 particles by default + central mass...
	assert len(sim.particles)==101
	ms=np.array([pp.m for pp in sim.particles[1:]])
	assert np.all(ms==5e-5)
	odat=sim.calculate_orbits(primary=sim.particles[0])
	eccs=np.array([oo.e for oo in odat])
	aa=np.array([oo.a for oo in odat])
	assert np.allclose(eccs, 0.7)
	tt= kstest(aa, cum_dist)
	assert tt.pvalue>0.05

def test_init_circ():
	bc.bash_command('rm archive*')
	bc.bash_command('python ../end_aleksey_config_b.py --config config2')
	name=bc.bash_command("echo archive*bin").replace('\n', '')
	sa=rebound.SimulationArchive(name)
	sim=sa[0]
	##100 particles by default + central mass...
	assert len(sim.particles)==201
	ms=np.array([pp.m for pp in sim.particles[1:]])
	assert np.all(ms==5e-4)
	odat=sim.calculate_orbits(primary=sim.particles[0])
	eccs=np.array([oo.e for oo in odat])
	assert np.allclose(eccs, 0.)


def test_frac():
	bc.bash_command('rm archive*')
	bc.bash_command('python ../end_aleksey_config_b.py --config config_frac')
	name=bc.bash_command("echo archive*bin").replace('\n', '')
	sa=rebound.SimulationArchive(name)
	sim=sa[0]
	ms=np.array([pp.m for pp in sim.particles[1:]])
	assert np.allclose(ms, 0.01/200.0)

def test_delete():
	bc.bash_command('rm archive*')
	bc.bash_command('python ../end_aleksey_config_b.py --config config_delete')
	name=bc.bash_command("echo archive*bin").replace('\n', '')
	sa=rebound.SimulationArchive(name)
	sim=sa[0]

	ms=np.array([pp.m for pp in sim.particles[1:]])
	assert len(ms[ms==1.0e-5])==100
	assert len(ms[ms==1.0e-4])==10


