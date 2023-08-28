# reproduce functionality of the currentTRDNR and powerTRDNR functions from Andreas' Mathematica
# notebook in Python

# https://doi.org/10.1103/PhysRevApplied.12.064018
# photon flux with mu and Eg:
from solcore.constants import kb, c, h, q
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *

def planck_dist(Eph, mu, kT):
    return Eph**2/(np.exp((Eph-mu)/kT)-1)

def planck_heaviside_dist(Eph, Eg, mu, kT):
    hvs_weight = np.heaviside(Eph-Eg, 0.5)
    pd_ys = planck_dist(Eph, mu, kT)
    return pd_ys*hvs_weight

def spec_pflux_planckheaviside(Ephs, Eg, mu, kT):
    return ((2 * np.pi) / (c ** 2 * (h / q) ** 3)) * planck_heaviside_dist(Ephs, Eg, mu, kT)

def Ndot_boltzmann(Eg, T, mu): # particle flux density
    # accurate for large band gaps / large negative bias
    kT = kb*T/q # convert to eV to match units of Eg
    N = ((2*np.pi)/(c**2*h**3))*np.exp((mu-Eg)/kT)*kT*(Eg**2 + 2*Eg*kT + 2*kT**2)
    # paper above, eq. 13
    return N*q**3


def Ndot_planckheaviside(Ephs, Eg, mu, kT):
    spec_pflux = spec_pflux_planckheaviside(Ephs, Eg, mu, kT)
    pflux = np.trapz(spec_pflux, Ephs)
    return pflux

def Ndot_downwellheaviside(Eg, downwell_dict):
    Ephs = downwell_dict['photon energies']
    get_heavisided = np.heaviside(Ephs - Eg, 0.5) * downwell_dict['downwelling photon flux']
    pflux = np.trapz(get_heavisided, Ephs)
    return pflux




Tc = 300  # temperature of emitter / converter
mu = 0
kT_c_eV = kb * Tc / q

Te = 30 # temperature of environment
kT_e_eV = kb * Te / q



Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
Egs = np.arange(0.01, 0.20, 0.002)

# fig, axs = plt.subplots(2,2, layout='tight')
#
# col_dict = {'planck_emit':'skyblue', 'planckheavy_emit':'darkblue', 'planck_env':'orangered', 'boltzmann_emit':'springgreen', 'downwell_env':'lightpink', 'downwellheavy_env':'deeppink'}
#
# Eg_single = 0.075  # [eV]
#
# spec_pflux_planck = ((2*np.pi)/(c**2*(h/q)**3))* planck_dist(Ephs, mu, kT_c_eV)  # [m-2.s-2/eV]
# spec_pflux_plankheaviside = spec_pflux_planckheaviside(Ephs, Eg_single, mu, kT_c_eV)  # [m-2.s-2/eV]
#
# axs[0][0].plot(Ephs, spec_pflux_planck, label='Planck dist, Blackbody at T=300K', color=col_dict['planck_emit'])
# axs[0][0].plot(Ephs, spec_pflux_plankheaviside, label=f'Planck dist $\\times$ Heaviside, $E_g$={Eg_single}eV', color=col_dict['planckheavy_emit'])
# axs[0][0].fill_between(Ephs, spec_pflux_plankheaviside, color=col_dict['planckheavy_emit'], alpha=0.1)
# axs[0][0].set_ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
# axs[0][0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
# axs[0][0].set_title('Photons emitted by converter per E$_{ph}$')
# axs[0][0].legend()
#
#
# N_boltzmann = Ndot_boltzmann(Egs, Tc, mu)
# N_planck_0Inf = np.vectorize(Ndot_planckheaviside, excluded=[0])(Ephs, Egs, mu, kT_c_eV)
#
# axs[0][1].plot(Egs, N_planck_0Inf, label='Planck $\\times$ Heaviside', color=col_dict['planckheavy_emit'])
# axs[0][1].plot(Egs, N_boltzmann, label='Boltzmann', color=col_dict['boltzmann_emit'])
# axs[0][1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
# axs[0][1].set_xlabel('Bandgap, E$_g$ [eV]')
# axs[0][1].set_title('Photons emitted by converter')
# axs[0][1].legend()
#
#
# # Environment
#
# # Spectral Photon Flux (per photon energy)
# pflux_env_blackbody = ((2*np.pi)/(c**2*(h/q)**3))* planck_dist(Ephs, mu, kT_e_eV)
# pflux_downwell = retrieve_downwelling_in_particleflux(10)
# Ephs_dw, spec_pflux_dw = pflux_downwell['photon energies'], pflux_downwell['downwelling photon flux']
# spec_pflux_dwheaviside = np.heaviside(Ephs_dw-Eg_single, 0.5)*spec_pflux_dw
#
# axs[1][0].plot(Ephs, pflux_env_blackbody, label=f'Blackbody, T={Te}K', color=col_dict['planck_env'])
# axs[1][0].plot(Ephs_dw, spec_pflux_dw, label='Downwelling (atmospheric modelling)', color=col_dict['downwell_env'])
# axs[1][0].plot(Ephs_dw, np.heaviside(Ephs_dw-Eg_single, 0.5)*spec_pflux_dw, label=f'Downwelling (absorbed for E$_g$ = {Eg_single} eV)', color=col_dict['downwellheavy_env'])
# axs[1][0].fill_between(Ephs_dw, spec_pflux_dwheaviside, color=col_dict['downwellheavy_env'], alpha=0.1)
#
# axs[1][0].set_ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
# axs[1][0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
# axs[1][0].set_title('Photons emitted by environment per E$_{ph}$')
# axs[1][0].legend()
#
#
# N_downwell_abs = np.vectorize(Ndot_downwellheaviside)(Egs, pflux_downwell)
#
# axs[1][1].plot(Egs, N_downwell_abs, label = 'Downwelling $\\times$ Heaviside', color=col_dict['downwellheavy_env'])
# axs[1][1].set_title('Photons abs from environment')
# axs[1][1].set_xlabel('Bandgap, E$_g$ [eV]')
# axs[1][1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
# axs[1][1].legend()
# plt.show()



# Drawing I-V, power curves
fig2, axs2 = plt.subplots(2,2, layout='tight')
cmap = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']
pflux_downwell = retrieve_downwelling_in_particleflux(10)
Egs = np.arange(0.062, 0.20, 0.002)

for iE, Eg in enumerate(Egs[::10]):
    # fermi level splitting
    mus = np.linspace(-Eg, 0, 100)

    # N_out, emitted from converter
    # N_in, emitted from environment, abs by converter
    N_out_mu_planck = np.vectorize(Ndot_planckheaviside, excluded=['Ephs'])(Ephs= Ephs, Eg=Eg, mu=mus, kT=kT_c_eV)
    N_in_mu_planck = Ndot_planckheaviside(Ephs, Eg, 0, kT_e_eV)
    N_in_mu_dw = Ndot_downwellheaviside(Eg, downwell_dict=pflux_downwell)

    # calculate and plot with N abs from env from (1/orange) blackbody at T=30K, (2/pink) downwelling data
    for i, N_in, cmap in zip([0,1],[N_in_mu_planck, N_in_mu_dw], [cmap_orange, cmap_pink]):
        i_colour = cmap(0.2 + 0.8 * iE / len(Egs[::10]))
        J_mu = q*(N_out_mu_planck - N_in)  # [J] [m-2.s-1] --> A.m-2
        # J = N_out - N_in

        axs2[i][0].plot(mus, J_mu, label = f'E$_g$ = {Eg:.3f} eV', color=i_colour)  # current
        axs2[i][1].plot(mus, mus*J_mu, color=i_colour)  # power density


for ri in [0,1]:
    for ci in [0,1]:
        axs2[ri][ci].plot([-1,1], [0,0], '-k', lw=0.8)
        axs2[ri][ci].set_xlim([-Egs[-1], 0])


    axs2[ri][0].set_xlabel("Fermi level splitting, $\mu$ (eV)")
    axs2[ri][0].set_ylabel("Current (A.m$^{-2}$)")
    axs2[ri][0].set_title('Current vs $\mu$')
    axs2[ri][0].legend()

    axs2[ri][1].set_xlabel("Fermi level splitting, $\mu$ (eV)")
    axs2[ri][1].set_ylabel("Power Density (W.m$^{-2}$)")
    axs2[ri][1].set_title('Power Density vs $\mu$')

plt.show()



# Limits: what changes with real atmospheric/device data?
# - non-radiative recombination (eq. 22)
# - N_environment isn't from a black body at a fixed temperature. Not all photons which are incident are absorbed,
# not all internally-generated photons are emitted (EQE/EL correspondence)
# - I think relevant quantity is the downwelling flux (possibly angle-dependent)? It doesn't matter if photons
# which are emitted by the device are absorbed later on in the atmosphere somewhere.


