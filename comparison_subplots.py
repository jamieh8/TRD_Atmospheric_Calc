import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

# cwv, T = 10, 296.724
# cwv, T = 24, 304.868
cwv, T = 54, 303.512

atm_data = atmospheric_dataset(cwv=cwv)
emitter_planck = planck_law_body(T=300)

Ephs = atm_data.photon_energies

angle_array = np.linspace(0,90,100)

comparison_lst = []

# comparing different cutoff angles

cmap = plt.get_cmap('tab10')
cutoff_angles = np.arange(10,95,20)
for ai, cutoff_angle in enumerate(cutoff_angles):
    comparison_lst += [{'label':f'cwv{cwv}, cutoff {cutoff_angle}', 'color':cmap(ai),
                        'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
                        'cutoff angle':cutoff_angle}]


# integrate over solid angle, using heaviside weighing and different angles of acceptance
# cutoff angle

include_photflux_plots = True
if include_photflux_plots:
    fig_pf, axs_pf = plt.subplots(2,2, layout='tight')

Egs = np.arange(0.062, 0.20, 0.002)

include_muopt_plots = True
if include_muopt_plots:
    fig_Popt, axs_Popt = plt.subplots(1,2, layout='tight')

for sample_dct in comparison_lst:
    atm_data = sample_dct['atmospheric dataset']
    emitter = sample_dct['emitter body']
    cutoff_angle = sample_dct['cutoff angle']

    if include_photflux_plots:
        # outgoing photon flux
        spec_phot_flux_out = emitter.spectral_photon_flux(Ephs, mu=0, cutoff_angle=cutoff_angle)
        axs_pf[0][0].plot(Ephs, spec_phot_flux_out, label=f'cutoff {cutoff_angle}', color=sample_dct['color'])

        Ndot_out = np.vectorize(emitter.retrieve_Ndot_heaviside, excluded=[0])(Ephs, Egs, mu=0, cutoff_angle=cutoff_angle)
        axs_pf[0][1].plot(Egs, Ndot_out, color=sample_dct['color'])

        # incoming photon flux
        spectral_photon_flux = atm_data.spectral_data_with_cutoffangle(angle_array, cutoff_angle)
        axs_pf[1][0].plot(Ephs, spectral_photon_flux, label=sample_dct['label'], color=sample_dct['color'])

        Ndot_vs_Eg = np.vectorize(atm_data.retrieve_Ndot_heaviside)(Egs, cutoff_angle)
        axs_pf[1][1].plot(Egs, Ndot_vs_Eg, label=sample_dct['label'], color=sample_dct['color'])


    if include_muopt_plots:
        # TRD performance metrics
        combined_obj = TRD_in_atmosphere(emitter, atm_data, Ephs)
        maxPs = []
        Vmpps = []
        for Eg in Egs:
            mu_opt = combined_obj.optimize_mu(Eg, cutoff_angle)
            maxPs += [mu_opt['max power']]
            Vmpps += [mu_opt['Vmpp']]
        axs_Popt[0].plot(Egs, maxPs, color=sample_dct['color'], label=sample_dct['label'])
        axs_Popt[1].plot(Egs, Vmpps, color=sample_dct['color'])

if include_photflux_plots:
    axs_pf[0][0].legend()
    axs_pf[1][0].legend()

    for spec_ax in [axs_pf[0][0], axs_pf[1][0]]:
        spec_ax.set_ylabel('Spectral Phot. D.F. [s$^{-1}$.m$^{-2}$/eV]')
        spec_ax.set_xlabel('Photon Energy, E$_{ph}$ [eV]')

    for Ndot_ax in [axs_pf[0][1], axs_pf[1][1]]:
        Ndot_ax.set_ylabel('Photon Density Flux, $\mathrm{\dot{N} \;[s^{-1}.m^{-2}]}$')
        Ndot_ax.set_xlabel('Bandgap, E$_g$ [eV]')

    axs_pf[0][1].set_title('Photons emitted by TRD')
    axs_pf[1][1].set_title('Photons absorbed from environment')

if include_muopt_plots:
    axs_Popt[0].set_xlabel('Bandgap, E$_g$ [eV]')
    axs_Popt[0].set_ylabel('Max Power Density [W.m$^{-2}$]')
    axs_Popt[0].legend()

    axs_Popt[1].set_xlabel('Bandgap, E$_g$ [eV]')
    axs_Popt[1].set_ylabel('V$_{mpp}$ [V]')


plt.show()
