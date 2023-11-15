import numpy as np
import refl_coef as rc
import liq_cryst as lc
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

# free space wavelength
lambda_0 = 1550 # in nanometers

# index of refraction
n_LC = (1.5 + 1.7)/2
n_H = 2.32
n_L = 1.45
n_SiO2 = 1.5
# vector with refraction indices
n = [n_LC, n_H, n_L, n_H, n_L, n_H, n_L, n_H, n_L, n_H, n_SiO2]

# primary reflection indices
rho_1 = (n_LC - n_H)/(n_LC + n_H)
rho_10 = (n_H - n_SiO2)/(n_H + n_SiO2)
rho_H2L = (n_H - n_L)/(n_H + n_L) # high to low 
rho_L2H = -rho_H2L                # low to high
# vector with primary reflection indices
rho = [rho_1, rho_H2L, rho_L2H, rho_H2L, rho_L2H, rho_H2L, rho_L2H, rho_H2L, rho_L2H, rho_10]

# width of layer with the higher index of refraction
d_H = lambda_0/4/n_H
# width of layer with the lower index of refraction
d_L = lambda_0/4/n_L 
# vector with the width of the layers
d = [d_H, d_L, d_H, d_L, d_H, d_L, d_H, d_L, d_H]

# knowledge of vectors rho, n, and d enable the calculation 
# of the reflection coefficient at any interface
lambdas = np.linspace(1000, 2000, 501)
refl_coeff = []
power_rc = []
abs_refl_coeff = []
for lambda_i in lambdas:
    r = rc.refl_coef(n, rho, d, lambda_i)
    refl_coeff.append(r)
    abs_refl_coeff.append(np.abs(r))
    power_rc.append(np.abs(r)**2)

plt.figure()
plt.plot(lambdas, 20*np.log10(abs_refl_coeff))
plt.ylabel("Refl. Coefficient(dB)")
plt.xlabel("Wavelength(nm)")
plt.xlim(1000,2000)
plt.grid()

plt.figure()
plt.plot(lambdas, 10*np.log10(power_rc))
plt.ylabel("Power Refl. Coefficient(dB)")
plt.xlabel("Wavelength(nm)")
plt.xlim(1000,2000)
plt.grid()

# Due to reciprocity, the power transmission coefficient
# is the same for either side of the Bragg reflector
# Since R+T = 1 AND R'+T' = 1 with T = T' => R = R'
# So both the power reflection coeff and the transmission coeff
# stay the same.
power_transm_c = 1 - np.array(power_rc)
power_rc = np.array(power_rc)  
refl_coeff = np.array(refl_coeff)
# The Liquid Crystal Layer has the following width
d_LC = 2*lambda_0/n_LC
# Normal incidence implies theta_i = 0 => cos(theta_i) = 1
phase_delta = np.divide(4 * np.pi * n_LC * d_LC, lambdas)

# Using Airy's relations with R = power_rc and T = 1 - power_rc
t_l = np.divide(np.abs(power_transm_c), np.abs(1-np.multiply(refl_coeff**2, np.exp(-1j*phase_delta))))
plt.figure()
plt.plot(lambdas, 20*np.log10(t_l))
plt.grid()
plt.title(r'$T(\lambda)$')
plt.ylabel(r'power transm. coeff. $10\log(|T(\lambda)|)$ (dB)')
plt.xlabel(r'Wavelength $\lambda$(nm)')
# Now instead of using Airy's formulae, we will use the path integral method
# taking into account every possible path an incident ray may follow
# As depicted in the image the incident wave comes from the SiO2 substrate and penetrates
# the multi-layer structure

# the d vector becomes
d = [*d, d_LC, *d]
# the n vector becomes
n_bragg = [n_H, n_L, n_H, n_L, n_H, n_L, n_H, n_L, n_H]
n = [n_SiO2, *n_bragg, n_LC, *n_bragg, n_SiO2]
# the rho vector becomes
rho = [-rho_10, rho_H2L, rho_L2H, rho_H2L, rho_L2H, rho_H2L, rho_L2H,
       rho_H2L, rho_L2H, -rho_1, *rho]

# calculate refl. coeff.
power_rc_alternative = []
for lambda_i in lambdas:
    r = rc.refl_coef(n, rho, d, lambda_i)
    power_rc_alternative.append(np.abs(r)**2)
power_transm_c_alternative = 1 - np.array(power_rc_alternative)

plt.plot(lambdas, 10*np.log10(power_transm_c_alternative), '*')
plt.legend(['A(ii)', 'A(iii)'])

######################################################################
# widths for the Bragg layers
d_bragg = [d_H, d_L, d_H, d_L, d_H, d_L, d_H, d_L, d_H]
rho_bragg = [rho_H2L, rho_L2H, rho_H2L, rho_L2H, rho_H2L, rho_L2H, 
             rho_H2L, rho_L2H]

# Liquid crystal layer width
d_LC = 1900 # in nanometers 

metadata = dict(title='convergence')
writer = FFMpegWriter(fps=1, metadata=metadata)

# The liquid crystal layer is divided into L sublayers
L = 14
lambdas = np.linspace(1450, 1700, 1001)
fig = plt.figure()
plt.xlim(1450,1700)

with writer.saving(fig, "fitPlot.mp4", 100):
    for L in range(2,15):    
        plt.clf()
        plt.grid(True)
        # for V = 1.5
        d_lc_1, n_lc_1, rho_lc_1 = lc.sublayer_properties(L, 1.5, 1.7, n_H, 1.5, d_LC)
        d1 = [*d_bragg, *d_lc_1, *d_bragg]
        n1 = [n_SiO2, *n_bragg, *n_lc_1, *n_bragg, n_SiO2]
        rho1 = [-rho_10, *rho_bragg, *rho_lc_1, *rho_bragg, rho_10]
        # for V = 3.0
        d_lc_2, n_lc_2, rho_lc_2 = lc.sublayer_properties(L, 1.5, 1.7, n_H, 3, d_LC)
        d2 = [*d_bragg, *d_lc_2, *d_bragg]
        n2 = [n_SiO2, *n_bragg, *n_lc_2, *n_bragg, n_SiO2]
        rho2 = [-rho_10, *rho_bragg, *rho_lc_2, *rho_bragg, rho_10]
        # for V = 4.5
        d_lc_3, n_lc_3, rho_lc_3 = lc.sublayer_properties(L, 1.5, 1.7, n_H, 4.5, d_LC)
        d3 = [*d_bragg, *d_lc_3, *d_bragg]
        n3 = [n_SiO2, *n_bragg, *n_lc_3, *n_bragg, n_SiO2]
        rho3 = [-rho_10, *rho_bragg, *rho_lc_3, *rho_bragg, rho_10]

        T1 = []
        T2 = []
        T3 = []
        for lambda_i in lambdas:
            r1 = rc.refl_coef(n1, rho1, d1, lambda_i) 
            r2 = rc.refl_coef(n2, rho2, d2, lambda_i) 
            r3 = rc.refl_coef(n3, rho3, d3, lambda_i) 
            T1.append(rc.power_transm_coeff(r1))
            T2.append(rc.power_transm_coeff(r2))
            T3.append(rc.power_transm_coeff(r3))
        T1_dB = 10*np.log10(T1)
        T2_dB = 10*np.log10(T2)
        T3_dB = 10*np.log10(T3)
        plt.plot(lambdas, T1_dB)
        plt.plot(lambdas, T2_dB)
        plt.plot(lambdas, T3_dB)
        plt.title(f'L = {L}')
        plt.xlabel(r'Wavelength $\lambda$ (nm)')
        plt.ylabel(r'Power Transm. Coeff. $10\log|T(\lambda)|$ (dB)')
        plt.legend(['V = 1.5 Volt','V = 3 Volt','V = 4.5 Volt'])


        # wavelength at resonance
        mx_ind1 = np.argmax(T1_dB)
        lambda_res_1 = lambdas[mx_ind1]
        # resonance BW
        tmp_array = lambdas[T1_dB >= -3]
        delta_lambda_3dB_1 = tmp_array[-1] - tmp_array[0]
        tmp_array = lambdas[T1_dB >= 10*np.log10(0.9)]
        delta_lambda_10dB_1 = tmp_array[-1] - tmp_array[0]

        # wavelength at resonance
        mx_ind2 = np.argmax(T2_dB)
        lambda_res_2 = lambdas[mx_ind2]
        # resonance BW
        tmp_array = lambdas[T2_dB >= -3]
        delta_lambda_3dB_2 = tmp_array[-1] - tmp_array[0]
        tmp_array = lambdas[T2_dB >= 10*np.log10(0.9)]
        delta_lambda_10dB_2 = tmp_array[-1] - tmp_array[0]

        # wavelength at resonance
        mx_ind3 = np.argmax(T3_dB)
        lambda_res_3 = lambdas[mx_ind3]
        # resonance BW
        tmp_array = lambdas[T3_dB >= -3]
        delta_lambda_3dB_3 = tmp_array[-1] - tmp_array[0]
        tmp_array = lambdas[T3_dB >= 10*np.log10(0.9)]
        delta_lambda_10dB_3 = tmp_array[-1] - tmp_array[0]

        print([lambda_res_1, delta_lambda_3dB_1, delta_lambda_10dB_1])
        print([lambda_res_2, delta_lambda_3dB_2, delta_lambda_10dB_2])
        print([lambda_res_3, delta_lambda_3dB_3, delta_lambda_10dB_3])
        
        writer.grab_frame()

# reveals all previous plots
plt.show()