import numpy as np
import math

#----two triangles--------
#parameters
mu=-1.0
kT=1.0
k_insertion=0.1

#list to save number-of-subunit values
Ns = []

#initial number of subunits
n=1
moves = ['insertion', 'removal']

#time steps loop
for t in xrange(2000000):
    #choose a move type randomly
    move = np.random.choice(moves)  
    if move=='insertion':
        if n==1:
            # 3.0x due to the number of edges
            p_propose = 3.0*k_insertion
            r_insert = np.random.rand()
            if r_insert<p_propose:
                p_acc = math.exp(mu/kT)
                r = np.random.rand()
                if r<p_acc:
                    n+=1

    if move=='removal':
        if n==2:
            # 2.0x because there are two subunits
            p_propose = 2.0*k_insertion
            r_remove = np.random.rand()
            if r_remove<p_propose:
                p_acc = math.exp(-mu/kT)
                r=np.random.rand()
                if r<p_acc:
                    n-=1
    Ns.append(n)
    
#compute free energy and entropy
N = np.arange(1,4)
P_N = np.histogram(Ns, bins=N)
P_N = P_N[0] / float(sum(P_N[0]))

#grand canonical potential
Omega_N = -kT*np.log(P_N)
Omega_N = Omega_N - Omega_N[0]

#Free energy
F_N = Omega_N+mu*N[:-1]
F_N = F_N - F_N[0]

#entropy
S_N = -F_N/kT
print "S(N)=", S_N

#number of states
W_N = np.exp(S_N)
print "W(N)=",W_N