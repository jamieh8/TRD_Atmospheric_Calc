# reproduce functionality of the currentTRDNR and powerTRDNR functions from Andreas' Mathematica
# notebook in Python

# https://doi.org/10.1103/PhysRevApplied.12.064018
# photon flux with mu and Eg:
from solcore.constants import kb, c, h, q
import numpy as np
import matplotlib.pyplot as plt

def planck_dist(Eph, mu, kT):
    return Eph**2/(np.exp((Eph-mu)/kT)-1)

def planck_heaviside_dist(Eph, Eg, mu, kT):
    hvs_weight = np.heaviside(Eph-Eg, 0.5)
    pd_ys = planck_dist(Eph, mu, kT)
    return pd_ys*hvs_weight

def flux_boltzmann(Eg, T, mu): # particle flux density
    # accurate for large band gaps / large negative bias
    kT = kb*T/q # convert to eV to match units of Eg
    N = ((2*np.pi)/(c**2*h**3))*np.exp((mu-Eg)/kT)*kT*(Eg**2 + 2*Eg*kT + 2*kT**2)
    # paper above, eq. 13
    return N*q**3

def flux_planck(Eg, T, mu, dE=0.0001, int_from_Eg=False):
    kT = kb*T/q
    if int_from_Eg:
        Eph = np.arange(Eg, 10*kT, dE)
        pd_ys = planck_dist(Eph, mu, kT)
    else:
        Eph = np.arange(1e-6, 10*kT, dE)
        pd_ys = planck_heaviside_dist(Eph, Eg, mu, kT)

    N = ((2*np.pi)/(c**2*h**3))*np.trapz(pd_ys, Eph)
    return N*q**3


T = 300  # temperature of environemt
mu = 0
kT = kb * T / q

fig, axs = plt.subplots(1,2, layout='tight')

Eg_single = 0.05
Ephs = np.arange(1e-6, 10*kT, 0.0001)
pd = planck_dist(Ephs, mu, kT)
pd_hv = planck_heaviside_dist(Ephs, Eg_single, mu, kT)
axs[0].plot(Ephs, pd, label='Planck dist, Blackbody at T=300K')
axs[0].plot(Ephs, pd_hv, '--', dashes=(3, 3), label=f'Planck dist $\\times$ Heaviside, $E_g$={Eg_single}eV')
axs[0].set_ylabel('Photon Density Flux, per energy [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
axs[0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
axs[0].set_title('Photons emitted per E$_{ph}$')
axs[0].legend()


Egs = np.linspace(0.01, 0.10, 50)
N_boltzmann = flux_boltzmann(Egs, T, mu)
N_planck_EgInf = np.vectorize(flux_planck, cache=False)(Egs, T, mu, int_from_Eg=True)
N_planck_0Inf = np.vectorize(flux_planck)(Egs, T, mu, int_from_Eg=False)
E = np.linspace(0, 1, 100)

axs[1].plot(Egs, N_planck_EgInf, label='Planck (int E$_g$ to $\\infty$)')
axs[1].plot(Egs, N_planck_0Inf, '--', dashes=(3, 3), label='Planck $\\times$ Heaviside (int 0 to $\\infty$)')
axs[1].plot(Egs, N_boltzmann, label='Boltzmann')
axs[1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
axs[1].set_xlabel('Bandgap, E$_g$ [eV]')
axs[1].set_title('Photons emitted for given E$_g$')
axs[1].legend()



# Drawing I-V, power curves
Eg = 1.5*kT
Te = 0.1*T  # temperature of environment, 30 K.

# fermi level splitting
mus = np.linspace(-Eg, 0, 100)

# N_out, emitted from converter
# N_in, emitted from environment, abs by converter
N_out_mu_planck = np.vectorize(flux_planck)(Eg, T, mus, int_from_Eg=False)
N_in_mu_planck = flux_planck(Eg, Te, 0, int_from_Eg=True)
J_mu_planck = q*(N_out_mu_planck - N_in_mu_planck)  # [J] [m-2.s-1] --> A.m-2



fig2, axs2 = plt.subplots(1,2, layout='tight')
axs2[0].plot(mus, J_mu_planck, label = 'Planck (emit & abs)')
axs2[0].set_xlabel("Fermi level splitting, $\mu$ (eV)")
axs2[0].set_ylabel("Current (A.m$^{-2}$)")
axs2[0].set_title('Current vs $\mu$')
axs2[0].legend()

# power density = J*V
axs2[1].plot(mus, mus*J_mu_planck)
axs2[1].set_xlabel("Fermi level splitting, $\mu$ (eV)")
axs2[1].set_ylabel("Power Density (W.m$^{-2}$)")
plt.show()



# Limits: what changes with real atmospheric/device data?
# - non-radiative recombination (eq. 22)
# - N_environment isn't from a black body at a fixed temperature. Not all photons which are incident are absorbed,
# not all internally-generated photons are emitted (EQE/EL correspondence)
# - I think relevant quantity is the downwelling flux (possibly angle-dependent)? It doesn't matter if photons
# which are emitted by the device are absorbed later on in the atmosphere somewhere.


