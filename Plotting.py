# reproduce functionality of the currentTRDNR and powerTRDNR functions from Andreas' Mathematica
# notebook in Python

# https://doi.org/10.1103/PhysRevApplied.12.064018
# photon flux with mu and Eg:
from solcore.constants import kb, c, h, q
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar


Tc = 300  # temperature of emitter / converter
mu = 0
kT_c_eV = kb * Tc / q

Te = 30 # temperature of environment
kT_e_eV = kb * Te / q
k_eV = kb/q



Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

pflux_downwell = retrieve_downwelling_in_particleflux(10)

col_dict = {'planck_emit':'skyblue', 'planckheavy_emit':'darkblue', 'planck_env_30K':'orangered', 'planck_env_100K':'orange',
            'boltzmann_emit':'springgreen', 'downwell_env':'lightpink', 'downwellheavy_env':'deeppink'}



fig, axs = plt.subplots(2,2, layout='tight')

Egs = np.arange(0.01, 0.20, 0.002)
Eg_single = 0.075  # [eV]

spec_pflux_planck = ((2*np.pi)/(c**2*(h/q)**3))* planck_dist(Ephs, mu, kT_c_eV)  # [m-2.s-2/eV]
spec_pflux_plankheaviside = spec_pflux_planckheaviside(Ephs, Eg_single, mu, kT_c_eV)  # [m-2.s-2/eV]

axs[0][0].plot(Ephs, spec_pflux_planck, label='Planck dist, Blackbody at T=300K', color=col_dict['planck_emit'])
axs[0][0].plot(Ephs, spec_pflux_plankheaviside, label=f'Planck dist $\\times$ Heaviside, $E_g$={Eg_single}eV', color=col_dict['planckheavy_emit'])
axs[0][0].fill_between(Ephs, spec_pflux_plankheaviside, color=col_dict['planckheavy_emit'], alpha=0.1)
axs[0][0].set_ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
axs[0][0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
axs[0][0].set_title('Photons emitted by converter per E$_{ph}$')
axs[0][0].legend()


N_boltzmann = Ndot_boltzmann(Egs, Tc, mu)
N_planck_0Inf = np.vectorize(Ndot_planckheaviside, excluded=[0])(Ephs, Egs, mu, kT_c_eV)

axs[0][1].plot(Egs, N_planck_0Inf, label='Planck $\\times$ Heaviside', color=col_dict['planckheavy_emit'])
axs[0][1].plot(Egs, N_boltzmann, label='Boltzmann', color=col_dict['boltzmann_emit'])
axs[0][1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
axs[0][1].set_xlabel('Bandgap, E$_g$ [eV]')
axs[0][1].set_title('Photons emitted by converter')
axs[0][1].legend()


# Environment

# Spectral Photon Flux (per photon energy)
pflux_env_30K = ((2*np.pi)/(c**2*(h/q)**3))* planck_dist(Ephs, mu, k_eV*30)
pflux_env_100K = ((2*np.pi)/(c**2*(h/q)**3))* planck_dist(Ephs, mu, k_eV*100)
Ephs_dw, spec_pflux_dw = pflux_downwell['photon energies'], pflux_downwell['downwelling photon flux']
spec_pflux_dwheaviside = np.heaviside(Ephs_dw-Eg_single, 0.5)*spec_pflux_dw

axs[1][0].plot(Ephs, pflux_env_30K, label=f'Blackbody, T=30K', color=col_dict['planck_env_30K'])
axs[1][0].plot(Ephs, pflux_env_100K, label=f'Blackbody, T=100K', color=col_dict['planck_env_100K'])
axs[1][0].plot(Ephs_dw, spec_pflux_dw, label='Downwelling (atmospheric modelling)', color=col_dict['downwell_env'])
axs[1][0].plot(Ephs_dw, np.heaviside(Ephs_dw-Eg_single, 0.5)*spec_pflux_dw, label=f'Downwelling (absorbed for E$_g$ = {Eg_single} eV)', color=col_dict['downwellheavy_env'])
axs[1][0].fill_between(Ephs_dw, spec_pflux_dwheaviside, color=col_dict['downwellheavy_env'], alpha=0.1)

axs[1][0].set_ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
axs[1][0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
axs[1][0].set_title('Photons emitted by environment per E$_{ph}$')
axs[1][0].legend()


# Photon Density Flux (Ndot per bandgap)
N_downwell_abs = np.vectorize(Ndot_downwellheaviside)(Egs, pflux_downwell)
N_planck30K = np.vectorize(Ndot_planckheaviside, excluded=[0])(Ephs, Egs, 0, k_eV*30)
N_planck100K = np.vectorize(Ndot_planckheaviside, excluded=[0])(Ephs, Egs, 0, k_eV*100)

axs[1][1].plot(Egs, N_downwell_abs, label = 'Downwelling $\\times$ Heaviside', color=col_dict['downwellheavy_env'])
axs[1][1].plot(Egs, N_planck30K, label = 'T=30K $\\times$ Heaviside', color=col_dict['planck_env_30K'])
axs[1][1].plot(Egs, N_planck100K, label = 'T=100K $\\times$ Heaviside', color=col_dict['planck_env_100K'])
axs[1][1].set_title('Photons abs from environment')
axs[1][1].set_xlabel('Bandgap, E$_g$ [eV]')
axs[1][1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
axs[1][1].legend()
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


    axs2[ri][0].set_xlabel("Fermi level splitting, $\mu$ [eV]")
    axs2[ri][0].set_ylabel("Current [A.m$^{-2}$]")
    axs2[ri][0].set_title('Current vs $\mu$')
    axs2[ri][0].legend()

    axs2[ri][1].set_xlabel("Fermi level splitting, $\mu$ [eV]")
    axs2[ri][1].set_ylabel("Power Density [W.m$^{-2}$]")
    axs2[ri][1].set_title('Power Density vs $\mu$')
#
# plt.show()


# Plotting max power density vs Eg
plt.figure()
# neg_powerdensity_downwellheaviside(mu, Eg, Ephs, kT_c_eV, pflux_downwell)
# find best mu for given Eg
Egs_dw = np.arange(0.062, 0.20, 0.002)
Egs_planck = np.arange(0, 0.20, 0.002)
maxPs_dwh = []
maxPs_30K = []
maxPs_100K = []

for Eg in Egs_dw:
    opt_mu_dwh = minimize_scalar(neg_powerdensity_downwellheaviside, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, pflux_downwell))
    maxPs_dwh += [opt_mu_dwh.fun]

for Eg in Egs_planck:
    opt_mu_30K = minimize_scalar(neg_powerdensity_plancks, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, k_eV*30))
    maxPs_30K += [opt_mu_30K.fun]

    opt_mu_100K = minimize_scalar(neg_powerdensity_plancks, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, k_eV * 100))
    maxPs_100K += [opt_mu_100K.fun]

plt.plot(Egs_dw, maxPs_dwh, color=col_dict['downwellheavy_env'], label='300K out, Atmospheric+Heaviside in')
plt.plot(Egs_planck, maxPs_30K, color=col_dict['planck_env_30K'], label='300K out, 30K in')
plt.plot(Egs_planck, maxPs_100K, color='orange', label='300K out, 100K in')
plt.ylabel('Power Density [W.m$^{-2}$]')
plt.xlabel('Bandgap [eV]')
plt.title('Max Power Density vs Bandgap')
plt.legend()
plt.show()


# Limits: what changes with real atmospheric/device data?
# - non-radiative recombination (eq. 22)
# - N_environment isn't from a black body at a fixed temperature. Not all photons which are incident are absorbed,
# not all internally-generated photons are emitted (EQE/EL correspondence)
# - I think relevant quantity is the downwelling flux (possibly angle-dependent)? It doesn't matter if photons
# which are emitted by the device are absorbed later on in the atmosphere somewhere.


