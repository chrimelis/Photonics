import numpy as np

# The functions below stem from the standard micro-ring coupling analysis
# supposing known: 1) losses per circulation a
#                  2) t1,t2 : from the lumped model that describes the coupling of the MRR
#                             and the access waveguide
#                  3) phi: phase accumulated per circulation (and more generally also due to t1,t2)

def T_through(t1, t2, a, phi):
    return np.abs((t2**2 * a**2 - 2*a*t1*t2*np.cos(phi) + t1**2) / (1 - 2*t1*t2*a*np.cos(phi) + (t1*t2*a)**2))

def T_drop(t1, t2, a, phi):
    return np.abs((1-t2**2)*(1-t1**2)*a / (1 - 2*t1*t2*a*np.cos(phi) + (t1*t2*a)**2))