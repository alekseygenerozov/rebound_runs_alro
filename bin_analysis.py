import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations


def bin_props(p1, p2):
	'''
	Auxiliary function to get binary properties for two particles. 

	p1 and p2 -- Two particles from a rebound simulation.
	'''
	dp = p1 - p2   
	##Calculate com velocity of two particles...
	##Masses
	m1=p1.m
	m2=p2.m
	com = (m1*p1+m2*p2)/(m1+m2)
	##Particle velocities in com frame
	p1_com = p1 - com
	p2_com = p2 - com
	v12 = (p1_com.vx**2.+p1_com.vy**2.+p1_com.vz**2.)
	v22 = (p2_com.vx**2.+p2_com.vy**2.+p2_com.vz**2.)

	d12 = (p1_com.x**2.+p1_com.y**2.+p1_com.z**2.)
	d22 = (p2_com.vx**2.+p2_com.y**2.+p2_com.z**2.)			

	##Kinetic and potential energies
	ke = 0.5*m1*v12+0.5*m2*v22
	d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
	##Assumes G has been set to 1.
	pe = (m1*m2)/d2**0.5

	##Distance of binary center of mass from COM of system (should be near central SMBH)
	com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
	a_bin=(m1*m2)/(2.*(pe-ke))
	j_bin=m1*np.cross(p1_com.xyz, p1_com.vxyz)+m2*np.cross(p2_com.xyz, p2_com.vxyz)
	mu=m1*m2/(m1+m2)

	e_bin=(1.-np.linalg.norm(j_bin)**2./((m1+m2)*a_bin)/(mu**2.))

	return com_d, a_bin, e_bin, p1_com, p2_com, d2

def bin_find(loc):
	##Ensure we are in the com frame of the simulation.
	t,name=loc
	sat = rebound.SimulationArchive(name)
	sim = sat.getSimulation(t)

	sim.move_to_com()
	ps = sim.particles
	##mass of of primary 
	m0 = ps[0].m
	bin_indics=[]
	for i1, i2 in combinations(range(1, sim.N),2): # get all pairs of indices/start at 1 to exclude SMBH
		com_d, a_bin, e_bin, p1_com, p2_com, d2 = bin_props(ps[i1], ps[i2])
		m1,m2 =ps[i1].m, ps[i2].m
		##Hill sphere condition.
		inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)

		##If the kinetic energy is less than the potential energy 
		if ((a_bin>0) and (inside_hill)):
			bin_indics.append([sim.t, i1, i2, d2**0.5, a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), e_bin])

	return bin_indics

def p_dist(loc, idx):
	t,name=loc
	sat = rebound.SimulationArchive(name)
	sim = sat.getSimulation(t)
	ps=sim.particles
	ds=np.empty(len(ps))
	for ii in range(0,len(ps)):
		dp=ps[ii]-ps[idx]
		ds[ii] = (dp.x*dp.x+dp.y*dp.y+dp.z*dp.z)**0.5

	order=np.argsort(ds)
	return order, ds[order]

def com_plot(sa_name, i1, i2, name=''):
	fig,ax=plt.subplots(figsize=(10,9))
	
	sa = rebound.SimulationArchive(sa_name)
	m0=sa[0].particles[0].m
	for ii,sim in enumerate(sa):
		p1,p2=sim.particles[i1],sim.particles[i2]
		m1,m2=p1.m,p2.m
		com_d, a_bin, e_bin, p1_com, p2_com, d2 = bin_props(p1,p2)
		p1_pos=p1_com.x, p1_com.y
		p2_pos=p2_com.x, p2_com.y

		ax.set_xlim(-0.1,  0.1)
		ax.set_ylim(-0.1, 0.1)

		ann=ax.annotate('a={0:2.2f}, a/rt={1:2.2f}, e^2={2:2.2f}'.format(a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), e_bin), (0.09, 0.09), horizontalalignment='right',\
							verticalalignment=20, fontsize=20)
		ax.plot(p1_pos[0], p1_pos[1], 'rs')
		ax.plot(p2_pos[0], p2_pos[1], 'ks')
		# print name+'com2_{0:03d}.png'.format(ii)
		fig.savefig(name+'com_{0:03d}.png'.format(ii))
		ann.remove()

class BinAnalysis(object):
	def __init__(self, sa_name):
		'''
		Getting properties of all of the binaries in a rebound simulation run.
		'''
		self.sa_name=sa_name
		sa = rebound.SimulationArchive(sa_name)

		self.m0=sa[0].particles[0].m
		self.ts= [sim.t for sim in sa]
		self.delta_t=np.diff(self.ts)[0]
		self.locs = [[sim.t, sa_name] for sim in sa]
		try:
			self.bins=np.genfromtxt(sa_name.replace('.bin','_bins.csv'), delimiter=',')
		except:
			self.__bin_init__()

		self.pairs_arr=self.bins[:,[1,2]].astype(int)
		self.times_arr=self.bins[:,0]
		self.pairs=np.array([{int(self.bins[i,1]),int(self.bins[i,2])} for i in range(len(self.bins))])



	def __bin_init__(self):
		pool = rebound.InterruptiblePool(processes=3)
		bins = pool.map(bin_find,self.locs)
		bins=np.array(bins)
		filt=[bins[i]!=[] for i in range(len(bins))]
		bins2=np.array(bins[filt])
		bins2=np.concatenate(bins2)
		self.bins=bins2

		np.savetxt(self.sa_name.replace('.bin','_bins.csv'), self.bins,delimiter=',')


	def exotica(self):
		'''
		Identifying triples and exchange interactions.
		'''
		for ns in range(np.min(self.pairs_arr), np.max(self.pairs_arr)+1):
			last=[]
			tlast=0.
			for tt in np.unique(self.times_arr):
				filt=(self.times_arr==tt)
				filt2=[(ns in pp) for pp in self.pairs[filt]]
				tmp=self.pairs[filt][filt2]

				if len(tmp)>1:
					print "star {0}, {1} bound stars!, tt={2}".format(ns, len(tmp)+1, tt/(2.*np.pi))
				elif (len(tmp)==1) and (len(last)==1) and (list(tmp)!=list(last)) and (tt-tlast<1.01*self.delta_t):
					print "star {0}, exchange!, tt={1}, {2}->{3}, {4}".format(ns, tt/(2.*np.pi), last, tmp, (tt-tlast)/self.delta_t)
				last=np.copy(tmp)
				tlast=tt

	def num_bins(self):
		'''
		Number of binaries for each snapshot of the simulation.
		'''
		num_bins=[len(self.times_arr[self.times_arr==tt]) for tt in self.ts]
		return num_bins






