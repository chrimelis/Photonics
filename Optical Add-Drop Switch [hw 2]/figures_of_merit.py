import numpy as np
import transm_coeff

def ER_max(T_port):
    ''' Input:  @ T_port : the transmission coeff regarding a specific port
                           for various frequencies/wavelengths/phase constants
                           
        Output: @ ER_max_dB: the maximum available extinction ratio 
                             for the given T_port matrix 
    '''
    if min(T_port) == 0:
        ER_max_dB = -10**4 
    else:
        ER_max_dB = 10*np.log10(min(T_port)/max(T_port)) 
    return ER_max_dB
    
def IL_min(T_port):
    ''' Input:  @ T_port : the transmission coeff regarding a specific port
                           for various frequencies/wavelengths/phase constants
                           
        Output: @ IL_min_dB: the minimum available insertion loss 
                             for the given T_port matrix 
    '''
    return 10*np.log10(1/max(T_port))

def merit_figures(phi,t1,t2,a):
    ''' Goal: For the given set of values (t1,t2,a) find Tthrough and Tdrop
        Input:  @ phi: = beta * 2pi R
                @ t1,t2: parameters of the coupling matrices
                @ a: losses per circulation
    '''
    T_through = [0]*len(phi)
    T_drop = [0]*len(phi)
    for k, phi_k in enumerate(phi):
        T_through[k] = transm_coeff.T_through(t1, t2, a, phi_k)
        T_drop[k] = transm_coeff.T_drop(t1, t2, a, phi_k)

    # Evaluation of the optical switch's performance
    # Figures of merit used: ER_max & IL_min at the output ports(drop & through)

    # Through port
    ER_through = ER_max(T_through)
    IL_through = IL_min(T_through) 
    # Drop port
    ER_drop = ER_max(T_drop)
    IL_drop = IL_min(T_drop)
    return (ER_through,ER_drop,IL_through,IL_drop)