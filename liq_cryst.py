import numpy as np

def theta_LC(z, V, d_LC):
    theta_max = 10 + 40/3 * V  #linear relation
    theta_max = theta_max * np.pi/180 # convert from deg to rad 
    return theta_max * np.sin(np.pi * z / d_LC) 

def n_eff(no, ne, z, V, d_LC):
    theta = theta_LC(z, V, d_LC)
    return no*ne/np.sqrt(no**2 * np.cos(theta)**2 + ne**2 * np.sin(theta)**2) 

def sublayer_properties(L, no, ne, n_H, V, d_LC):
    ''' Inputs:
    @ L: the number of sublayers
    @ no: the ordinary refr. index
    @ ne: the extrordinary refr. index
    @ V: the voltage applied
    @ d_LC: the width of the liquid crystal layer
        Outputs:
    1) the width of each sublayer (which is equal to d_LC / L) as list with length L
    2) the refractive index vector
    3) the principal reflection coefficient for each sublayer
    '''
    # 1)
    d_sub = d_LC / L
    # 2) 
    n = []
    for i in range(L):
        # For the sake of simplicity, the value of neff is determined by the
        # height of the middle plane of the sublayer
        # of the sublayer. 
        # In the limit as L->oo the choice of the height is irrelevant
        z_tmp = d_sub/2 + i * d_sub
        n_tmp = n_eff(no, ne, z_tmp, V, d_LC)
        n.append(n_tmp)
    # 3)
    rho = [ (n_H - n[0]) / (n_H + n[0]) ]
    for i in range(L-1):
        rho_tmp = (n[i] - n[i+1])/(n[i] + n[i+1])
        rho.append(rho_tmp)
    rho_tmp = (n[-1] - n_H)/(n[-1] + n_H)
    rho.append(rho_tmp)
    return [d_sub]*L, n, rho