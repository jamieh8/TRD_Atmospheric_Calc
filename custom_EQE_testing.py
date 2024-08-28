import numpy as np
from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

# load output from raytracing as 2 arrays (angles, escape_prob(angles))
# escape_prob = np.loadtxt('sphere_raytrace_totalT_2e15_70_points_4900_rays.txt', skiprows=0)
# betas_mod = np.radians(np.linspace(0,90,100))


# load measured EQE from csv as 2 arrays (wavelengths, EQE(wavelengths))
# fn = r"C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Temp dependent EQE\EQE_5um_diode.csv"
fn = r'C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Andreas Dropbox - Calculations\EQE_5um_diode.csv'
data = np.loadtxt(fn, skiprows=3, delimiter=',')
wavelengths = data[:,0]
EQE_PVI5 = data[:,3]/100

# load ray tracing results as 2 arrays (modelled emission angles beta, escape_prob(beta)
fn = r"C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Andreas Dropbox - Calculations\hyperhemi_T_withinangle_u_45.txt"
raytrace_data = np.loadtxt(fn, dtype=float)
betas_mod = raytrace_data[:,0]
escape_prob = raytrace_data[:,1]

# set up the diode using measured EQE and simulated escape probability
PVI5 = diode_EQE(spec_type='measured', spec_xs={'spec x':wavelengths, 'spec units':'wavelength [um]'}, spec_EQE=EQE_PVI5,
                            angle_type='escape probability - isotropic', escape_prob=escape_prob, betas=betas_mod,
                            resistance=377, area=0.1*1e-6)  # area 0.1mm2


# What is the power output of this diode at 350K, if the environment is 300K ?

# set up object that defines the emission spectrum / photon flux...
m350 = mathemetica_style_body(T=350, nL=3.3, R_normal=1-0.714)

# emitted photon flux from this diode, at short circuit:
Ndot_emitted_SC = m350.retrieve_Ndot(PVI5, mu=0)  # photon flux density emitted at short circuit (mu=qV=0)

# assuming absorption from the 300K environment
m300 = mathemetica_style_body(T=300, nL=3.3, R_normal=1-0.714)
Ndot_absorbed = m300.retrieve_Ndot(PVI5, mu=0)

# calculating current from photon flux
Jsc = 2 * q * (Ndot_emitted_SC - Ndot_absorbed)  # factor of 2: emission and absorption across full sphere, not just hemisphere
print(Jsc)  # result 581.3 (Python) compared to 580.973 (Mathematica)


# alternatively, use the diode_in_environment class which combines emitter and environment
diode5_Tc350_Te300 = diode_in_environment(emitter=m350, environment=m300, diode_EQE=PVI5)
Jsc = 2 * diode5_Tc350_Te300.current_density(V=0)  # factor of 2 added to consider a full sphere
print(Jsc)  # same result as above, but easier to call

# plotting an IV
# Vs = np.arange(0, -0.1, -0.005)
# print(Vs)
# Js, PDs = [], []
# for V in Vs:
#     Js += [2 * diode5_Tc350_Te300.current_density(V=V)]
#     PDs += [2 * diode5_Tc350_Te300.power_density(V=V)]
# plt.plot([0,-0.1], [0,0], '-k')
# plt.xlim([-0.1,0])
# # plt.plot(Vs, Js)
# plt.plot(Vs, PDs)

# checking power outputs
m3 = mathemetica_style_body(T=3, nL=3.3, R_normal=1-0.714)
m293 = mathemetica_style_body(T=293, nL=3.3, R_normal=1-0.714)
diode5_Tc293_Te3 = diode_in_environment(emitter=m293, environment=m3, diode_EQE=PVI5)
print(diode5_Tc293_Te3.powerdensity_from_R())

# check env temp vs power density
# Ts = np.arange(100, 292, 1)
# maxPDs, maxPDs_withPb = [], []
# for env_T in Ts:
#     menv = mathemetica_style_body(T=env_T)
#     die = diode_in_environment(emitter=m293, environment=menv, diode_EQE=PVI5)
#     mpp = die.powerdensity_from_R()
#     maxPDs += [mpp]
#
# plt.plot(Ts, 1e3*np.array(maxPDs), label='PVI5, hemisphere')
# # plt.plot(Ts, np.array(maxPDs)/min(maxPDs))
# plt.xlabel('Temp of environment [K]')
# plt.ylabel('Power density [mW.m-2]')
# plt.ylabel('Power density (relative)')


plt.show()