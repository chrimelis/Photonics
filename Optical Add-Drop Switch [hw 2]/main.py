import numpy as np
import matplotlib.pyplot as plt
import pi_ticks as pticks
import transm_coeff
import figures_of_merit as fig_of_m
import optimizers as opt

# Add-Drop Filter implemented with 1 ring resonator and
# 2 access waveguides as depicted below
###########################
# Drawing of our device:
#   ----------
#       O
#   ----------
###########################

#---------------------------------------------------------
# Part 1: Performance of the switch at critical coupling
#---------------------------------------------------------
# Parameters of the coupling Matrices
t1 = 0.9
t2 = 0.9
# losses per circulation
a = 0.85
# exponential decay of the form exp(-ar * theta * R)
# Hence, exp(-ar * 2pi * R) = a
# <=>    ar = -ln(a)/(2*pi*R)
#
# For 0.5 circulations the loss factor will be:
# exp(-ar * pi * R) = exp(lna/(2*pi*R) * pi * R) 
# = sqrt(a)

# phi = beta * 2*pi*R
phi = np.linspace(np.pi, 5*np.pi, 501)

# Transmission coefficients:
    # T_through = (t2**2 * a**2 - 2*a*t1*t2*cos(phi) + t1**2) / (1 - 2*t1*t2*a*cos(phi) + (t1*t2*a)**2)
    # T_drop = (1-t2**2)*(1-t1**2)*a / (1 - 2*t1*t2*a*cos(phi) + (t1*t2*a)**2)
T_through = [0]*len(phi)
T_drop = [0]*len(phi)

for k, phi_k in enumerate(phi):
    T_through[k] = transm_coeff.T_through(t1, t2, a, phi_k)
    T_drop[k] = transm_coeff.T_drop(t1, t2, a, phi_k)

T_through_dB = 10*np.log10(T_through)
T_drop_dB = 10*np.log10(T_drop)

plt.figure()
plt.plot(phi, T_through_dB)
plt.xlabel(r'$\phi = \beta 2 \pi  R$ (radians)')
plt.ylabel(r'power transm. coeff. $10\log(T_{through})$ (dB)')
plt.title('Through Transmission Coeff.')
plt.grid(True)
ax = plt.gca()
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.xaxis.set_major_formatter(plt.FuncFormatter(pticks.multiple_formatter()))

plt.figure()
plt.plot(phi, T_drop_dB)
plt.xlabel(r'$\phi = \beta 2 \pi  R$ (radians)')
plt.ylabel(r'power transm. coeff. $10\log(T_{drop})$ (dB)')
plt.title('Drop Transmission Coeff.')
plt.grid(True)
ax = plt.gca()
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.xaxis.set_major_formatter(plt.FuncFormatter(pticks.multiple_formatter()))

plt.figure()
plt.plot(phi,T_through_dB, phi, T_drop_dB)
plt.xlabel(r'$\phi = \beta 2 \pi  R$ (radians)')
plt.ylabel(r'power transm. coeff. (dB)')
plt.title('Transmission Coefficients')
plt.grid(True)
plt.legend(['Through', 'Drop'])
ax = plt.gca()
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.xaxis.set_major_formatter(plt.FuncFormatter(pticks.multiple_formatter()))

# Extinction Ratio Calculations
# For the through port:
ER_through = min(T_through_dB) - max(T_through_dB)
# For the drop port:
ER_drop = min(T_drop_dB) - max(T_drop_dB)
print(f"ER_through = {ER_through} dB & ER_drop = {ER_drop} dB.")

# Insertion Loss Calculations
# For the through port:
IL_through = 0 - max(T_through_dB)
# For the drop port:
IL_drop = 0 - max(T_drop_dB)
print(f"IL_through = {IL_through} dB & IL_drop = {IL_drop} dB.")
print('---------------')


# -------------------------------------------------------------
#    Part 2: Optimizing the Design of the 1x2 optical switch    
# -------------------------------------------------------------
a = 0.85    # Constant parameter

# We have to set some criteria for our design
# By choosing values for t1,t2 these criteria
# must be satisfied. 

# -------------------------------
#       Brute-Force Search
# -------------------------------

t_values = np.linspace(0.05, 0.95, 100) 
grid_lines = len(t_values)

ER_through = np.zeros((grid_lines,grid_lines))
IL_through = np.zeros((grid_lines,grid_lines))
ER_drop = np.zeros((grid_lines,grid_lines))
IL_drop = np.zeros((grid_lines,grid_lines))

for m, t2 in enumerate(t_values):
    for n, t1 in enumerate(t_values):
        # First for the given pair (t1,t2) find Tthrough and Tdrop
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


t1,t2 = np.meshgrid(t_values,t_values)

plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(t1,t2,ER_through)
plt.title('ER through(dB)')
plt.xlabel('t1')
plt.ylabel('t2')

plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(t1,t2,ER_drop)
plt.title('ER drop(dB)')
plt.xlabel('t1')
plt.ylabel('t2')

plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(t1,t2,IL_through)
plt.grid(True)
plt.title('IL through(dB)')
plt.xlabel('t1')
plt.ylabel('t2')

plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(t1,t2,IL_drop)
plt.grid(True)
plt.title('IL drop(dB)')
plt.xlabel('t1')
plt.ylabel('t2')

#------------------------------------------------------------------
#   No pairs (t1,t2) can satisfy simultaneously all target values
#------------------------------------------------------------------
condition_A = ER_through <= -15
condition_B = ER_drop <= -14
condition_C = IL_through <= 3
condition_D = IL_drop <= 5
# Combine conditions using logical AND
combined_condition = np.logical_and.reduce([condition_A, condition_B, condition_C, condition_D])
# Find the indices where all conditions are satisfied
indices = np.argwhere(combined_condition)
print('---------------')
print(f"Attempt for satisfaction of all criteria.\n indices = {indices}")  # => EMPTY SET
print('---------------')

#---------------------------------------------------------------
#   Weighted Cost Function for approximately acceptable design
#---------------------------------------------------------------

# Find the indices with the minimum cost
min_cost_indices = opt.cost_function(ER_through,ER_drop,IL_through,IL_drop)
t1 = t_values[min_cost_indices[1]]
t2 = t_values[min_cost_indices[0]]
ER_through, ER_drop, IL_through, IL_drop = fig_of_m.merit_figures(phi,t1,t2,a)
# print figures of merit
print('---------------')
print(f"(t1,t2) = ({t1,t2})")
print(f"ER_through = {ER_through} dB & ER_drop = {ER_drop} dB.")
print(f"IL_through = {IL_through} dB & IL_drop = {IL_drop} dB.")
print('---------------')

#--------------------------------------------------------------
#  Search in a strip around the critical coupling condition
#                    t1 = a*t2 +(or)- 0.1
#--------------------------------------------------------------
# The idea behind this search is that ER_through is known to 
# behave ideally at exactly t1 = a*t2 and thus our criteria
# will be met in an area around this line in the parameter
# space (t1,t2) with 0<t1<1 and 0<t2<1
 
print('---------------')
opt.critical_strip_search(a,phi) #=> prints ours design parameters and figures of merit achieved
print('---------------')

plt.show()
