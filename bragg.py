import numpy as np

def bragg_reflector(d_H, d_L, n_H, n_L, M):
    '''Implements (HL)^N H
    and computes its properties'''
    # compute the principal reflection coefficients
    # M-1 are fully determined
    rho_H2L = (n_H - n_L)/(n_H + n_L)
    rho_L2H = -rho_H2L
    
    d_bragg = []
    n_bragg = []
    rho_bragg = []    
    
    N = round((M-1)/2)
    for i in range(N):
        d_bragg.append(d_H)
        d_bragg.append(d_L)
        n_bragg.append(n_H)
        n_bragg.append(n_L)
        rho_bragg.append(rho_H2L)
        rho_bragg.append(rho_L2H)
            
    d_bragg.append(d_H)
    n_bragg.append(n_H)
    
    return d_bragg, n_bragg, rho_bragg

