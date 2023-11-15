import numpy as np

def recursive_formula(n_1, rho_1, d_1, Gamma_2, lambda_0):
    k0 = 2*np.pi/lambda_0; 
    return ( rho_1 + Gamma_2 * np.exp(-1j*2*k0*n_1*d_1) ) / ( 1 + rho_1*Gamma_2*np.exp(-1j*2*k0*n_1*d_1) )

def refl_coef(n, rho, d, lambda_0):
    M = len(d)
    r = rho[M]
    for i in range(M):
        r = recursive_formula(n[M-i], rho[M-1-i], d[i], r, lambda_0)
    return r

def power_transm_coeff(r):
    '''@ r: the reflection coefficient refering 
    to the ratio of reflected to incident electric field'''
    return 1-np.abs(r)**2