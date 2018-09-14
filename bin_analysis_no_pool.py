import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations


def get_com(ps):
	'''
	Get center of mass for a collection of particles
	'''
	com=ps[0].m*ps[0]
	ms=ps[0].m
	for pp in ps[1:]:
		com=com+pp.m*pp
		ms+=pp.m
	return com/ms

def bin_props(p1, p2):
	'''
	Auxiliary function to get binary properties for two particles. 

	p1 and p2 -- Two particles from a rebound simulation.
	'''
	dp = p1 - p2   
	d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
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
	##Difference in the forces acting on the two particles;
	ft = np.array([m2*(p2.ax)-m2*(com.ax), m2*(p2.ay)-m2*(com.ay), m2*(p2.az)-m2*com.az])
	##Unit vector pointing from particle 2 to particle 1
	rhat = np.array(dp.xyz)/d2**0.5
	f12 = m1*m2/d2*rhat 
	##Subtract the mutual force between the two stars
	ft = ft - f12
	ft = np.linalg.norm(ft)

	##Kinetic and potential energies
	ke = 0.5*m1*v12+0.5*m2*v22
	##Assumes G has been set to 1.
	pe = (m1*m2)/d2**0.5

	##Distance of binary center of mass from COM of system (should be near central SMBH)
	com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
	a_bin=(m1*m2)/(2.*(pe-ke))
	j_bin=m1*np.cross(p1_com.xyz, p1_com.vxyz)+m2*np.cross(p2_com.xyz, p2_com.vxyz)
	j_com=(m1+m2)*np.cross(com.xyz, com.vxyz)

	inc=np.arccos(np.dot(j_bin, j_com)/np.linalg.norm(j_bin)/np.linalg.norm(j_com))*180./np.pi
	mu=m1*m2/(m1+m2)
	e_bin=(1.-np.linalg.norm(j_bin)**2./((m1+m2)*a_bin)/(mu**2.))

	return com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft


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
		com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft = bin_props(ps[i1], ps[i2])
		m1,m2 =ps[i1].m, ps[i2].m
		##Hill sphere condition.
		inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)
		tidal_2 = (m1*m2/d2>ft)

		##If the kinetic energy is less than the potential energy 
		if ((a_bin>0) and (inside_hill) and (tidal_2)):
		#if ((a_bin>0) and (inside_hill)):
			bin_indics.append([sim.t, i1, i2, d2**0.5, a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), e_bin])

	return np.array(bin_indics)


def bin_find_sim(sim):
	##Ensure we are in the com frame of the simulation.
	sim.move_to_com()
	ps = sim.particles
	##mass of of primary 
	m0 = ps[0].m
	bin_indics=[]
	for i1, i2 in combinations(range(1, sim.N),2): # get all pairs of indices/start at 1 to exclude SMBH
		com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft = bin_props(ps[i1], ps[i2])
		m1,m2 =ps[i1].m, ps[i2].m
		##Hill sphere condition.
		inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)
		tidal_2 = (m1*m2/d2>ft)

		##If the kinetic energy is less than the potential energy 
		if ((a_bin>0) and (inside_hill) and (tidal_2)):
			bin_indics.append([sim.t, i1, i2, d2**0.5, a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), e_bin])

	return np.array(bin_indics)

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

def com_plot(sa_name, i1, i2, extras=[], name='', cols=['r', 'g', 'k'], idx_min=0, idx_max=None, lim=0.1):
	planes = [['x', 'y'], ['x', 'z'], ['y','z']]
	fig,ax=plt.subplots(nrows=1, ncols=3, figsize=(10*len(planes),10))
	#fig.patch.set_visible(False)

	for kk in range(3):
		ax[kk].set_xlim(-lim,  lim)
		ax[kk].set_ylim(-lim, lim)
		ax[kk].set_xlabel(planes[kk][0])
		ax[kk].set_ylabel(planes[kk][1])
		#ax[kk].patch.set_visible(False)

	sa = rebound.SimulationArchive(sa_name)
	if not idx_max:
		idx_max=len(sa)
	m0=sa[0].particles[0].m
	for ii in range(idx_min, idx_max):
		for kk, plane in enumerate(planes):
			sim=sa[ii]
			p1,p2=sim.particles[i1],sim.particles[i2]
			m1,m2=p1.m,p2.m
			com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft = bin_props(p1,p2)
			p1_pos=getattr(p1_com, plane[0]), getattr(p1_com, plane[1])
			p2_pos=getattr(p2_com, plane[0]), getattr(p2_com, plane[1])
			com=get_com([p1, p2])

			ax[kk].plot(p1_pos[0], p1_pos[1], 'o', markersize=2, color=cols[0])
			ax[kk].plot(p2_pos[0], p2_pos[1], 'o', markersize=2, color=cols[1])
			for jj, extra in enumerate(extras):
				ax[kk].plot(getattr(sim.particles[extra]-com, plane[0]), getattr(sim.particles[extra]-com, plane[1]), 'o', markersize=2, color=cols[(2+jj)%len(cols)])
			# if inset:
			# 	ax2.set_xlim(-(1+e_bin**0.5)*a_bin, (1+e_bin**0.5)*a_bin)
			# 	ax2.set_ylim(-(1+e_bin**0.5)*a_bin, (1+e_bin**0.5)*a_bin)
			# 	ax2.plot(p1_pos[0], p1_pos[1], 'o', markersize=2, color=cols[0])
			# 	ax2.plot(p2_pos[0], p2_pos[1], 'o', markersize=2, color=cols[1])

				# for jj,extra in enumerate(extras):
				# 	ax.plot(sim.particles[extra].x-com.x, sim.particles[extra].y-com.y, 'o', markersize=2, color=cols[(2+jj)%len(cols)])
			# print name+'com2_{0:03d}.png'.format(ii)
		ann=ax[1].annotate('a={0:2.2g}, a,/rt={1:2.2g}, r={2:2.2g}\n e^2={3:2.2g}, 1-e^2={4:2.2g}\n i={5:2.2g}'\
			.format(a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), com_d, e_bin, 1-e_bin, inc),\
			(0.9*lim, 0.9*lim), horizontalalignment='right',verticalalignment='top', fontsize=20)
		#fig.canvas.print_png(open(sa_name.replace('.bin', '')+name+'_com_{0:03d}.png'.format(ii), 'w'))
		fig.savefig(sa_name.replace('.bin', '')+name+'_com_{0:03d}.png'.format(ii), bbox_inches='tight', pad_inches=0)
		ann.remove()

class BinAnalysis(object):
	def __init__(self, sa_name):
		'''
		Getting properties of all of the binaries in a rebound simulation run.
		'''
		self.sa_name=sa_name
		sa = rebound.SimulationArchive(sa_name)

		self.m0=sa[0].particles[0].m
		self.ts= np.array([sim.t for sim in sa])
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
		#pool = rebound.InterruptiblePool(processes=3)
		bins = map(bin_find,self.locs[:10])
		#bins=np.array(bins)
		filt=np.array([len(bins[i])>0 for i in range(len(bins))])
		bins2=np.array(bins)[filt]
                print bins2
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

	def bin_times(self):
		pairs=self.pairs
		t_survs=np.zeros(len(np.unique(pairs)))
		sa = rebound.SimulationArchive(self.sa_name)
		##For each binary identify how long it survives
		for ii,pp in enumerate(np.unique(pairs)):
			##Identify all times where each binary pair exists.
			t_bin=self.times_arr[pairs==pp]
			t_surv=t_bin[-1]-t_bin[0]
			#Index of one of the stars in the pair
			idx=self.pairs_arr[pairs==pp][0,1]
			idx2=self.pairs_arr[pairs==pp][0,0]

			##Edge case: Binary splits up and forms again. See if the binary has skipped any snapshots.
			##snapshots may not be exactly evenly spaced, so it is best to to use the to if the binary 
			##has skipped any snapshots...
			diffs=np.diff([np.where(self.ts==t_bin[ii])[0][0] for ii in range(len(t_bin))])
			if np.any(diffs>1.01):
				tmp=np.split(t_bin, np.where(diffs>1.01)[0]+1)
				tmp2=[tmp[i][-1]-tmp[i][0] for i in range(len(tmp))]
				order=np.argsort(tmp2)
				t_bin=tmp[order[-1]]
				t_surv=tmp2[order[-1]]

			##Survival time of the binary normalized to orbital period (of binary within disk)
			##(we use the orbital period of one of the stars in the binary as a proxy). 
			sim=sa.getSimulation(t_bin[-1])
			sim.move_to_com()
			#t_orb=sim.particles[idx].P
			##Binary orbital period -- not this is not a constant--take the minimum orbital period 
			m1=sa[0].particles[idx].m
			m2=sa[0].particles[idx2].m
			t_orb = 2.*np.pi*np.min((self.bins[self.pairs==pp][:,4]**3./(m1+m2))**0.5)
			t_survs[ii]=t_surv/t_orb
		return t_survs







                
