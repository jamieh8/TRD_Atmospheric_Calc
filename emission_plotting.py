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

col_dict = {'planck_emit':'skyblue', 'planckheavy_emit':'darkblue', 'planck_env_30K':'orangered', 'planck_env_100K':'orange',
            'boltzmann_emit2':'springgreen', 'downwell_env':'lightpink', 'downwellheavy_env':'deeppink'}


# ------------ EFFECT OF MU ON TRD EMISSION CURVE ------------
fig, axs = plt.subplots(1,1, layout='tight', sharey='all')

Egs = np.arange(0.01, 0.20, 0.05)
Eg_single = 0.075  # [eV]

# Emissions from different mu
mus = [-0.075, -0.05, -0.025, 0]
ax = axs
mu_spec = -0.025
cmap = matplotlib.colormaps['plasma']
emitter = planck_law_body(T=300, Ephs = Ephs)
for mi, mu in enumerate(mus):
    spec_pflux_planck = emitter.spectral_photon_flux(Ephs, mu, cutoff_angle=90) #[s-1.m-2/eV]

    mu_colour = cmap(0.2 + 0.8 * mi / len(mus))
    ax.plot(Ephs, spec_pflux_planck, label=f'$\mu$ = {mu}', color=mu_colour)

    if mu == mu_spec:
        # spec_pflux_plankheaviside = emitter.retrieve_Ndot_heaviside(Eg=Eg_single, cutoff_angle=90, mu=mu_spec) #[s-1.m-2]
        # ax_area.plot(Ephs, spec_pflux_plankheaviside, label=f'Planck dist $\\times$ Heaviside, $E_g$={Eg_single}eV', color=mu_spec_colour)
        ax.fill_between(Ephs, spec_pflux_planck * np.heaviside(Ephs-Eg_single, 0.5), color=mu_colour, alpha=0.2)


ax.set_ylabel('SPFD, F$_\mathrm{ph}$ [$\mathrm{s^{-1}.m^{-2}.eV^{-1}}$]')
ax.set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
ax.set_title(f'Planck\'s law, T={Tc}K - Photons emitted by converter')
ax.legend()


# Plot downwelling data (vs photon energy)
# downwell_dat_dict = retrieve_downwelling_in_particleflux(10)
# plt.figure()
# plt.plot(downwell_dat_dict['photon energies'], downwell_dat_dict['downwelling photon flux'], c=col_dict['downwellheavy_env'])
# plt.plot([0.094,0.094], [0,4*1e23], '--k')
# plt.text(x=0.094, y=3*1e23, s='  0.094 eV',ha='left')
# plt.xlabel('Photon Energy, E$_{ph}$ [eV]')
# plt.ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}.eV^{-1}}$]')
#
# plt.ylim([0,3.8*1e23])



# Downwelling - W.m-2 vs Eph, to photon flux vs Eph
# atm_dat = atmospheric_dataset(10)
#
# fig, axs = plt.subplots(1,2, layout='tight')
# Ephs = atm_dat.photon_energies
#
# downwelling_pd = atm_dat.retrieve_spectral_array(yvals='W.m-2', xvals='eV', col_name='downwelling_flux')
# downwelling_pf = atm_dat.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV', col_name='downwelling_flux')
# axs[0].plot(Ephs, downwelling_pd, c='crimson')
# axs[1].plot(Ephs, downwelling_pf, c=col_dict['downwellheavy_env'])
#
# axs[0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
# axs[1].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
# axs[0].set_ylabel('(Spectral) Power Density Flux [$\mathrm{W.m^{-2}/eV}$]')
# axs[1].set_ylabel('(Spectral) Photon Density Flux [$\mathrm{m^{-2}.s^{-1}/eV}$]')


plt.show()