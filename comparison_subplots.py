import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

def Eph_to_sx(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def sx_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')


comparison_lst = []

# comparing different cutoff angles
# cwv, T = 10, 296.724
# cwv, T = 24, 304.868
# cwv, T = 54, 303.512
# atm_data = atmospheric_dataset(cwv=cwv)
# Ephs = atm_data.photon_energies
# emitter_planck = planck_law_body(T=T, Ephs=Ephs)
angle_array = np.linspace(0,90,100)
#
# cmap = plt.get_cmap('tab10')
# cutoff_angles = np.arange(10,95,20)
# for ai, cutoff_angle in enumerate(cutoff_angles):
#     comparison_lst += [{'label':f'cwv{cwv}, cutoff {cutoff_angle}', 'color':cmap(ai),
#                         'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
#                         'cutoff angle':cutoff_angle}]

# comparing different cwvs
Egs_AD = np.arange(0.062, 0.2, 0.002)
datsets = [{'cwv':10, 'Tc':296.724, 'colour':'deeppink'},
           {'cwv':24, 'Tc':304.868, 'colour':'steelblue'},
           {'cwv':54, 'Tc':303.512, 'colour':'seagreen'}]
for ds in datsets:
    cwv = ds['cwv']
    atm_data = atmospheric_dataset(cwv=cwv)
    Ephs = atm_data.photon_energies
    emitter_planck = planck_law_body(T=ds['Tc'], Ephs=Ephs)
    comparison_lst += [{'label': f'cwv{cwv}', 'color': ds['colour'],
                        'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
                         'cutoff angle':90, 'use diffusivity approx':True, 'Egs':Egs_AD}]


# comparing blackbody environments
Ephs = np.arange(1e-6, 0.31, 0.0001)
Egs_bb = np.arange(0,0.2,0.002)
emitter_planck_300 = planck_law_body(T=300, Ephs=Ephs)
Tsets = [{'Tc':3, 'colour':'orangered'}, {'Tc':200, 'colour':'orange'}]
for Ts in Tsets:
    Tc = Ts['Tc']
    bb_env = planck_law_body(Tc, Ephs)
    comparison_lst += [{'label': f'{Tc}K', 'color': Ts['colour'],
                        'atmospheric dataset': bb_env, 'emitter body': emitter_planck_300,
                        'cutoff angle': 90, 'use diffusivity approx': True, 'Egs':Egs_bb}]


include_photflux_plots = True
include_Ndot_diff = False
include_muopt_plots = False
include_heaviside_ex = False

if include_photflux_plots:
    fig_pf, axs_pf = plt.subplots(2,2, layout='tight')
    if include_Ndot_diff:
        fig_Ndotd, axs_Ndotd = plt.subplots(1,1, layout='tight')
if include_muopt_plots:
    fig_Popt, axs_Popt = plt.subplots(1,2, layout='tight')


for sample_dct in comparison_lst:
    atm_data = sample_dct['atmospheric dataset']
    emitter = sample_dct['emitter body']
    cutoff_angle = sample_dct['cutoff angle']
    Egs = sample_dct['Egs']

    if include_photflux_plots:
        # outgoing photon flux
        Ephs = emitter.Ephs
        spec_phot_flux_out = emitter.spectral_photon_flux(Ephs, mu=0, cutoff_angle=cutoff_angle)
        axs_pf[0][0].plot(Ephs, spec_phot_flux_out, label=sample_dct['label'], color=sample_dct['color'])

        Ndot_out = np.vectorize(emitter.retrieve_Ndot_heaviside)(Egs, cutoff_angle, 0)
        axs_pf[0][1].plot(Egs, Ndot_out, color=sample_dct['color'])

        # incoming photon flux
        if sample_dct['use diffusivity approx']:
            try:
                spectral_photon_flux = atm_data.retrieve_spectral_array(yvals = 's-1.m-2', xvals = 'eV', col_name = 'downwelling_flux')
            except:
                spectral_photon_flux = atm_data.spectral_photon_flux(Ephs, mu=0, cutoff_angle=cutoff_angle)
        else:
            spectral_photon_flux = atm_data.spectral_data_with_cutoffangle(angle_array, cutoff_angle)
        axs_pf[1][0].plot(Ephs, spectral_photon_flux, label=sample_dct['label'], color=sample_dct['color'])

        Ndot_in = np.vectorize(atm_data.retrieve_Ndot_heaviside)(Egs, cutoff_angle)
        axs_pf[1][1].plot(Egs, Ndot_in, label=sample_dct['label'], color=sample_dct['color'])

        if include_Ndot_diff:
            axs_Ndotd.plot(Egs, Ndot_out-Ndot_in, label=sample_dct['label'], color=sample_dct['color'])

        if include_heaviside_ex:
            Eg_ex = 0.08
            # shade area under spectral PDF that would be integrated
            axs_pf[0][0].fill_between(Ephs, spec_phot_flux_out*np.heaviside(Ephs-Eg_ex,1), color=sample_dct['color'], alpha=0.2)
            axs_pf[1][0].fill_between(Ephs, spectral_photon_flux*np.heaviside(Ephs-Eg_ex,1), color=sample_dct['color'], alpha=0.2)

            # plot point in Ndot curve corresponding to integral
            axs_pf[0][1].plot([Eg_ex], [emitter.retrieve_Ndot_heaviside(Eg_ex, 0, 90)], 'o', color=sample_dct['color'],alpha=0.2)
            axs_pf[1][1].plot([Eg_ex], [atm_data.retrieve_Ndot_heaviside(Eg_ex,90)], 'o',color=sample_dct['color'], alpha=0.2)


    if include_muopt_plots:
        # TRD performance metrics
        combined_obj = TRD_in_atmosphere(emitter, atm_data)
        maxPs = []
        Vmpps = []
        for Eg in Egs:
            mu_opt = combined_obj.optimize_mu(Eg, cutoff_angle)
            maxPs += [mu_opt['max power']]
            Vmpps += [mu_opt['Vmpp']]
        axs_Popt[0].plot(Egs, maxPs, color=sample_dct['color'], label=sample_dct['label'])
        axs_Popt[1].plot(Egs, Vmpps, color=sample_dct['color'])



if include_photflux_plots:
    # add reference line
    for ri in [0,1]:
        for ci in [0,1]:
            ylim_max = axs_pf[ri][ci].get_ylim()[1]
            E_ref = 0.16
            axs_pf[ri][ci].plot(2*[E_ref],[0,1.1*ylim_max], '--k')
            axs_pf[ri][ci].text(x=E_ref+0.001, y=ylim_max/2, s=f'{E_ref} eV')
            axs_pf[ri][ci].set_ylim([0,ylim_max])

    # add heaviside demo
    Eg = 0.08  # [eV]
    first_plt_dct = comparison_lst[0]

    axs_pf[0][0].legend()
    axs_pf[1][0].legend()

    for spec_ax in [axs_pf[0][0], axs_pf[1][0]]:
        spec_ax.set_ylabel('Spectral P.F.D., F$_{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
        spec_ax.set_xlabel('Photon Energy, E$_{ph}$ [eV]')
        secax = spec_ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
        secax.set_xlabel('Wavenumber [cm$^{-1}$]')
        secax.tick_params(axis='x', length=7)

    for Ndot_ax in [axs_pf[0][1], axs_pf[1][1]]:
        Ndot_ax.set_ylabel('Photon Flux Density, $\mathrm{\dot{N} \;[s^{-1}.m^{-2}]}$')
        Ndot_ax.set_xlabel('Bandgap, E$_g$ [eV]')
        secax = Ndot_ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
        secax.set_xlabel('Wavenumber [cm$^{-1}$]')

    axs_pf[0][1].set_title('Photons emitted by TRD')
    axs_pf[1][1].set_title('Photons absorbed from environment')

if include_muopt_plots:
    axs_Popt[0].set_xlabel('Bandgap, E$_g$ [eV]')
    axs_Popt[0].set_ylabel('Max Power Density [W.m$^{-2}$]')
    axs_Popt[0].legend()

    axs_Popt[1].set_xlabel('Bandgap, E$_g$ [eV]')
    axs_Popt[1].set_ylabel('V$_{mpp}$ [V]')

    for ax in axs_Popt:
        secax = ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
        secax.set_xlabel('Wavenumber [cm$^{-1}$]')

plt.show()
