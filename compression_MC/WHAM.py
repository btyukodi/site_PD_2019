import numpy as np

def iterate_wall_WHAM(N_jN, A_j, U_jN):
    jaxis=0
    Naxis=1
    S1_N = np.sum(N_jN, axis=jaxis)
    nj = (np.sum(N_jN, axis=Naxis)).astype(float)

    S2_N = np.sum( (np.exp(A_j-U_jN.T)*nj).T, axis=jaxis)
    if 0 in S2_N:
        print "S2_N =0 !"

    exp_A_j_new = np.sum(S1_N/S2_N * np.exp(-U_jN), axis=Naxis)    
    exp_A_j_new[exp_A_j_new==0]=1
    A_j_new =-np.log( exp_A_j_new )
    #S2_N = np.sum( (np.exp(A_j-U_jN.T+U_jN)*nj).T, axis=jaxis)
    
    #A_j_new =-np.log( np.sum(S1_N/S2_N, axis=Naxis))    
    
    return A_j_new

def get_P_N(N_jN, A_j, U_jN):
    jaxis=0
    Naxis=1
    S1_N = np.sum(N_jN, axis=jaxis)
    nj = np.sum(N_jN, axis=Naxis).astype(float)

    S2_N = np.sum( (np.exp(A_j-U_jN.T)*nj).T, axis=jaxis) 
    return S1_N/S2_N

def get_free_energy_vs_N(shelf):
    N_jN, A_j, U_jN = get_wall_WHAM_arrays(shelf)
    for i in range(50000):
        A_j_prev = A_j.copy()
        A_j = iterate_wall_WHAM(N_jN, A_j, U_jN)
    P_N = get_P_N(N_jN, A_j, U_jN)
    #P_N[P_N==0]=1.0
    P_N = P_N[1:]
    print np.mean((A_j-A_j_prev)**2)
    return -np.log(P_N), A_j


def get_wall_WHAM_arrays(shelf):
    depth = 2000.0
    N_jN = np.zeros(len(shelf.itervalues().next()[0]))
    U_jN = np.zeros(len(shelf.itervalues().next()[0]))
    for k in np.sort(np.array(shelf.keys()).astype(int)):
        h = shelf[str(k)]
        #print k, sum(h[0]), sum(h[0]>0)
        if sum(h[0])>0:
            N_jN = np.vstack((N_jN, h[0]))
            U_this = np.zeros(len(h[0]))+depth
            U_this[k:k+4] = 0.0
            U_jN = np.vstack((U_jN, U_this))
        else:
            break
            
    #Ncut = np.max(np.where(np.sum(N_jN, axis=0)>0)[0])-5    
    #Ncut = np.max(np.where(np.sum(U_jN==0, axis=1)==4))
    #N_jN = N_jN[1:, :Ncut]
    #U_jN = U_jN[1:, :Ncut]
    
    
    N_jN = N_jN[1:]
    U_jN = U_jN[1:]    
    #!N_jN = N_jN[3:,4:]
    #!U_jN = U_jN[3:,4:]
    
    cut = np.max(np.where(N_jN[-1]>0)[0])+1
    N_jN = N_jN[:,:cut]
    U_jN = U_jN[:,:cut]
    
    A_j = (np.random.rand(len(N_jN))-0.5)
    #N_jN: P(N) histogram counts of the j-th simulation
    #U_jN: The U=U(N), the bias potential of the j-th simulation
    return N_jN, A_j, U_jN