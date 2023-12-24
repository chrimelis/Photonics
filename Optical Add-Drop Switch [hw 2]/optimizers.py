import numpy as np
import matplotlib.pyplot as plt
import transm_coeff
import figures_of_merit as fig_of_m

def critical_strip_search(a,phi):
    
    num1 = 50
    num2 = 50
    
    ER_through = np.zeros((num2,num1))
    IL_through = np.zeros((num2,num1))
    ER_drop = np.zeros((num2,num1))
    IL_drop = np.zeros((num2,num1))

    T_through = [0]*len(phi)
    T_drop = [0]*len(phi)
    
    t2_values = np.linspace(0.05,0.95,num2)
    for m,t2 in enumerate(t2_values):
        t1_min = np.maximum(a*t2 - 0.1,0)
        t1_max = np.minimum(a*t2 + 0.1,1)
        t1_values = np.linspace(t1_min,t1_max,num1)
        for n,t1 in enumerate(t1_values):
            for k, phi_k in enumerate(phi):
                T_through[k] = transm_coeff.T_through(t1, t2, a, phi_k)
                T_drop[k] = transm_coeff.T_drop(t1, t2, a, phi_k)

            # Evaluation of the optical switch's performance
            # Figures of merit used: ER_max & IL_min at the output ports(drop & through)

            # Through port
            ER_through[m,n] = fig_of_m.ER_max(T_through)
            IL_through[m,n] = fig_of_m.IL_min(T_through) 
            # Drop port
            ER_drop[m,n] = fig_of_m.ER_max(T_drop)
            IL_drop[m,n] = fig_of_m.IL_min(T_drop)
    
    min_cost_ind = cost_function(ER_through,ER_drop,IL_through,IL_drop)
    t1,t2 = strip_search_indices_to_t1_t2(a,num1,min_cost_ind[0],min_cost_ind[1],t2_values)
    for k, phi_k in enumerate(phi):
        T_through[k] = transm_coeff.T_through(t1, t2, a, phi_k)
        T_drop[k] = transm_coeff.T_drop(t1, t2, a, phi_k)

    # Evaluation of the optical switch's performance
    # Figures of merit used: ER_max & IL_min at the output ports(drop & through)

    # Through port
    ER_through = fig_of_m.ER_max(T_through)
    IL_through = fig_of_m.IL_min(T_through) 
    # Drop port
    ER_drop = fig_of_m.ER_max(T_drop)
    IL_drop = fig_of_m.IL_min(T_drop)
    # print figures of merit
    print(f"(t1,t2) = ({t1,t2})")
    print(f"ER_through = {ER_through} dB & ER_drop = {ER_drop} dB.")
    print(f"IL_through = {IL_through} dB & IL_drop = {IL_drop} dB.")
            
def cost_function(ER_through,ER_drop,IL_through,IL_drop):
    ''' Given the matrices ER_through, ER_drop, IL_through, IL_drop
        a penalty is heuristically calculated for each matrix element
        that does not satisfy the imposed requirements.
        Finally the indices which minimize the cost function are returned
    '''
    penalty1 = np.maximum((ER_through+15)/np.max(np.abs(ER_through)), 0);  w1 = 10
    penalty2 = np.maximum((ER_drop+14)/np.max(np.abs(ER_drop)), 0);        w2 = 10
    penalty3 = np.maximum((IL_through-3)/np.max(np.abs(IL_through)), 0);   w3 = 1
    penalty4 = np.maximum((IL_drop-5)/np.max(np.abs(IL_drop)), 0);         w4 = 1

    cost = w1*penalty1 + w2*penalty2 + w3*penalty3 + w4*penalty4
    # Find the indices with the minimum cost
    min_cost_indices = np.unravel_index(np.argmin(cost), cost.shape)
    return min_cost_indices

def strip_search_indices_to_t1_t2(a,num1,m,n,t2_values):
    ''' This method is used only for the critical_strip_search optimizer
    
        Inputs: @ a : losses per circulation of the ring resonator
                @ num1: num of columns
                @ num2: num of rows
                @ t2_values: the values we are examining in the search of the (sub)optimal solution
        
        Ouputs: @t1: parameter of the coupling matrix corresponding to col n
                @t2: parameter of the coupling matrix corresponding to row m
    '''
    t2 = t2_values[m]
    t1_min = np.maximum(a*t2 - 0.1,0)
    t1_max = np.minimum(a*t2 + 0.1,1)
    t1_values = np.linspace(t1_min,t1_max,num1)
    t1 = t1_values[n]
    return (t1,t2)
