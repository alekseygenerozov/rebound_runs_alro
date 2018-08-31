from itertools import combinations


def bin_find(sim):
    ##Ensure we are in the com frame of the simulation.
    sim.move_to_com()
    ps = sim.particles
    m0 = ps[0].m
    i=0
    bin_indics=[]
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
        
        com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
        a_bin=(sim.G*m1*m2)/(2.*(pe-ke))
        ##Hill sphere condition.
        inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)

        ##If the kinetic energy is less than the potential energy 
        if ((ke<pe) and (inside_hill)):
            i+=1
            bin_indics.append([i1, i2, d2**0.5, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), a_bin, inside_hill, ke<pe])
    return bin_indics

def td(sim, i1, i2):
    ##Ensure we are in the com frame of the simulation.
    sim.move_to_com()
    ps = sim.particles
    m0 = ps[0].m
    
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
    d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
    
    ke=0.5*m1*v12+0.5*m2*v22
    pe=sim.G*(m1*m2)/d2**0.5
    com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
    a_bin=(sim.G*m1*m2)/(2.*(pe-ke))
    
    print a_bin/(((m1+m2)/m0)**(1./3.)*com_d)
    inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)

    com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
    return [i1, i2,  a_bin/(((m1+m2)/m0)**(1./3.)*com_d), inside_hill, ke<pe]
   