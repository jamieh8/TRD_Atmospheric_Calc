import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def Eph_to_v(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def v_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')


comparison_lst = []

# comparing different cutoff angles
# # # cwv, T = 10, 296.724
# # # cwv, T = 24, 304.868
# # cwv, T = 54, 303.512
# cwv_str = 'low'
# loc_str = 'telfer'
# Tskin = 299.86
# atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str)
# Ephs = atm_data.photon_energies
# emitter_planck = planck_law_body(T=Tskin, Ephs=Ephs)
# angle_array = np.arange(0,91,1)
# Egs_AD = np.arange(0.02, 0.2, 0.05)
#
# cmap = plt.get_cmap('tab10')
# cutoff_angles = np.arange(10,95,20)
# for ai, cutoff_angle in enumerate(cutoff_angles):
#     line_format_dct = {'color':cmap(ai), 'linestyle':'solid'}
#     comparison_lst += [{'label':f'{loc_str} {cwv_str}, cutoff {cutoff_angle}', 'color':cmap(ai),
#                         'line format':line_format_dct, 'scatter format':{'marker':'o', 'c':cmap(ai)},
#                         'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
#                         'cutoff angle':cutoff_angle, 'Egs':Egs_AD, 'use diffusivity approx':False}]
#
# line_format_diff = {'color':'black', 'linestyle':'dashed', 'dashes':(4,4)}
# comparison_lst += [{'label': f'{loc_str} {cwv_str}, diffusivity approx 53$^\circ \\times \pi$',
#                     'line format':line_format_diff, 'scatter format':{'marker':'o', 'c':'k'},
#                     'atmospheric dataset': atm_data, 'emitter body': emitter_planck,
#                     'cutoff angle': None, 'Egs': Egs_AD, 'use diffusivity approx': True}]


# comparing new datasets:
angle_array = [0]
Egs_AD = np.arange(0.0125, 0.3, 0.002)
datasets = [
    {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o'},
    # {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o'},
    # {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o'},

    # {'loc':'california', 'cwvstring':'low', 'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'mid', 'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'high', 'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'symbol': 's'},
    #
    # {'loc':'tamanrasset', 'cwvstring':'low', 'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'mid', 'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'high', 'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'symbol':'^'}
    ]

for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str, Tskin=ds['Tskin'])
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['color'], 'linestyle': 'solid'}
    scatter_format = {'c': ds['color'], 'marker': ds['symbol'], 'markersize':8}
    emitter_planck = planck_law_body(T=ds['Tskin'], Ephs=Ephs)
    comparison_lst += [{'label': f'{loc_str.capitalize()} {cwv_str}', 'line format':line_format_dct, 'scatter format':scatter_format,
                        'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
                         'cutoff angle':None, 'use diffusivity approx':True, 'Egs':Egs_AD}]


# comparing different cwvs
# Egs_AD = np.arange(0.062, 0.3, 0.002)
# datsets = [{'cwv':10, 'Tc':296.724, 'colour':'deeppink'},
#            # {'cwv':24, 'Tc':304.868, 'colour':'steelblue'},
#            # {'cwv':54, 'Tc':303.512, 'colour':'seagreen'}
#            ]
# for ds in datsets:
#     cwv = ds['cwv']
#     atm_data = atmospheric_dataset(cwv=cwv)
#     Ephs = atm_data.photon_energies
#     line_format_dct = {'color': ds['colour'], 'linestyle': 'solid'}
#     emitter_planck = planck_law_body(T=ds['Tc'], Ephs=Ephs)
#     comparison_lst += [{'label': f'cwv{cwv}', 'line format':line_format_dct,
#                         'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
#                          'cutoff angle':None, 'use diffusivity approx':True, 'Egs':Egs_AD}]


# BBs at T skin
Tsets = []
for ds in datasets:
    Tsets += [{'Tc':ds['Tskin'], 'colour':ds['color']}]

# comparing blackbody environments
# Ephs = np.arange(1e-6, 0.31, 0.0001)
# Egs_bb = np.arange(0.001,0.3,0.002)
# emitter_planck_300 = planck_law_body(T=300, Ephs=Ephs)
# # Tsets = [{'Tc':3, 'colour':'black'}]
# Tsets = [{'Tc':200, 'colour':'navy'}, {'Tc':270, 'colour':'blueviolet'}, {'Tc':290, 'colour':'mediumorchid'}]
#
# # effective temperatures
# for dataset_entry in comparison_lst:
#     Teffective = dataset_entry['atmospheric dataset'].effective_skytemp(300)
#     Tsets += [{'Tc':Teffective, 'colour':dataset_entry['line format']['color']}]
#
#
# for Ts in Tsets:
#     Tc = Ts['Tc']
#     bb_env = planck_law_body(Tc, Ephs)
#     line_format_dct = {'color': Ts['colour'], 'linestyle': 'dashed'}
#     comparison_lst += [{'label': 'T$_\mathrm{atm}$ = '+f'{Tc:.5g}K', 'line format':line_format_dct,
#                         'atmospheric dataset': bb_env, 'emitter body': emitter_planck_300,
#                         'cutoff angle': None, 'use diffusivity approx': True, 'Egs':Egs_bb}]
#
# comparison_lst[-1].update({'label position':'below'})  # required for Tatm Telfer high, to accomodate Telfer mid

#
# custom_muopt_legend = [Line2D([0],[0], color = 'k', linestyle='solid', label='LBLTRM modelling'),
#                        Line2D([0],[0], color = 'k', linestyle='dashed', label='Effective temperature approx'),
#                        Patch(facecolor='darkorange', label='Telfer low'),
#                        Patch(facecolor='darkviolet', label='Telfer mid'),
#                        Patch(facecolor='teal', label='Telfer high'),
#                        Patch(facecolor='black', label='3K BB')]

custom_muopt_legend = []

log_atmdat = True

include_photflux_plots = True
include_Ndot_diff = False
include_heaviside_ex = True
Eg_ex = 0.25

include_muopt_plots = False
opt_Eg_and_mu = False  # adds points at optimal Eg
log_power = True
atmdat_background = True
use_cust_legend = False

include_Eg_PD_scatter = False


secondary_ticks = 'wavenumber'


alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

if include_photflux_plots:
    fig_pf, axs_pf = plt.subplots(2,2, layout='tight')
    if include_Ndot_diff:
        fig_Ndotd, axs_Ndotd = plt.subplots(1,1, layout='tight')

if include_muopt_plots:
    fig_Popt, axs_Popt = plt.subplots(1,2, layout='tight')

    if atmdat_background:
        ref_dataset = atmospheric_dataset_new('low', 'telfer', 301.56)
        downwelling_photflux = ref_dataset.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV',
                                                                   col_name='downwelling_flux')
        dwn_flux_yaxs = axs_Popt[0].twinx()
        dwn_flux_yaxs.plot(ref_dataset.photon_energies, downwelling_photflux, c='lightgrey')
        dwn_flux_yaxs.set_ylabel('Spectral Photon Flux Density, $\mathrm{F_{ph} \; [s^{-1}.m^{-2}/eV]}$', color='lightgrey')
        dwn_flux_yaxs.tick_params(axis='y', labelcolor='lightgrey')
        dwn_flux_yaxs.text(s='Telfer low', x=0.09, y=0.5*1e23, ha='right', color='lightgrey')
        if log_atmdat:
            dwn_flux_yaxs.set_yscale('log')

if include_Eg_PD_scatter:
    fig_scatter, axs_scatter = plt.subplots(1, 1, layout='tight')

for sample_dct in comparison_lst:
    atm_data = sample_dct['atmospheric dataset']
    emitter = sample_dct['emitter body']
    cutoff_angle = sample_dct['cutoff angle']
    Egs = sample_dct['Egs']
    style_args = sample_dct['line format']
    try:
        linestyle = sample_dct['linestyle']
    except:
        linestyle = '-'


    if include_photflux_plots:
        # outgoing photon flux
        Ephs = emitter.Ephs
        spec_phot_flux_out = emitter.spectral_photon_flux(Ephs, mu=0, cutoff_angle=cutoff_angle)
        axs_pf[0][0].plot(Ephs, spec_phot_flux_out, label=sample_dct['label'], **style_args)

        Ndot_out = np.vectorize(emitter.retrieve_Ndot_heaviside)(Egs, cutoff_angle, 0)
        axs_pf[0][1].plot(Egs, Ndot_out, **style_args)

        # incoming photon flux
        if sample_dct['use diffusivity approx']:
            try:
                spectral_photon_flux = atm_data.retrieve_spectral_array(yvals = 's-1.m-2', xvals = 'eV', col_name = 'downwelling_flux')
            except:
                spectral_photon_flux = atm_data.spectral_photon_flux(Ephs, mu=0, cutoff_angle=cutoff_angle)
        else:
            spectral_photon_flux = atm_data.spectral_PDF_with_cutoffangle(angle_array, cutoff_angle)
        axs_pf[1][0].plot(Ephs, spectral_photon_flux, label=sample_dct['label'], **style_args)

        Ndot_in = np.vectorize(atm_data.retrieve_Ndot_heaviside)(Egs, cutoff_angle)
        axs_pf[1][1].plot(Egs, Ndot_in, label=sample_dct['label'], **style_args)

        if include_Ndot_diff:
            axs_Ndotd.plot(Egs, Ndot_out-Ndot_in, label=sample_dct['label'], **style_args)

        if include_heaviside_ex:
            # shade area under spectral PDF that would be integrated
            axs_pf[0][0].fill_between(Ephs, spec_phot_flux_out*np.heaviside(Ephs-Eg_ex,1), **style_args, alpha=0.2)
            axs_pf[1][0].fill_between(Ephs, spectral_photon_flux*np.heaviside(Ephs-Eg_ex,1), **style_args, alpha=0.2)

            # plot point in Ndot curve corresponding to integral
            axs_pf[0][1].plot([Eg_ex], [emitter.retrieve_Ndot_heaviside(Eg=Eg_ex, cutoff_angle=90, mu=0)], 'o', **style_args,alpha=0.2)
            axs_pf[1][1].plot([Eg_ex], [atm_data.retrieve_Ndot_heaviside(Eg_ex,90)], 'o',**style_args, alpha=0.2)

    combined_obj = TRD_in_atmosphere(emitter, atm_data)
    if include_muopt_plots:
        # TRD performance metrics
        maxPs = []
        Vmpps = []
        for Eg in Egs:
            mu_opt = combined_obj.optimize_mu(Eg, cutoff_angle)
            maxPs += [mu_opt['max power']]
            Vmpps += [mu_opt['Vmpp']]

        # if log_power:
        maxPs = (-1)*np.array(maxPs)

        axs_Popt[0].plot(Egs, maxPs, **style_args, label=sample_dct['label'])
        axs_Popt[1].plot(Egs, Vmpps, **style_args)

    if opt_Eg_and_mu and (include_muopt_plots or include_Eg_PD_scatter):
        # optimize over Eg and mu simulaneously
        opt_xs, opt_pd = get_best_pd(combined_obj, args_to_opt=['Eg','mu'],
                                         args_to_fix={'cutoff_angle':None, 'consider_nonrad':False, 'eta_ext':1}, alg=alg_de)

        pd = opt_pd[0] * (-1)
        print(pd)
        if log_power:
           y_offset = 0.15*pd
        else:
            y_offset = 0.1

        if include_muopt_plots:
            axs_Popt[0].plot(opt_xs['Eg'], pd, 'o', **style_args)
            if use_cust_legend != True:
                # if no custom legend, label optimal pts
                label_pos = sample_dct.get('label position')
                if label_pos == 'below':
                    axs_Popt[0].text(s = sample_dct['label'], x=opt_xs['Eg'], y=pd-y_offset, c=style_args['color'], ha='right', va='top')
                else:
                    axs_Popt[0].text(s=sample_dct['label'], x=opt_xs['Eg'], y=pd + y_offset, c=style_args['color'],
                                         ha='right', va='bottom')

        if include_Eg_PD_scatter:
            print(sample_dct['label'])
            print(pd)
            axs_scatter.plot(opt_xs['Eg'], pd, **sample_dct['scatter format'])




if include_photflux_plots:
    # add reference line
    # for ri in [0,1]:
    #     for ci in [0,1]:
    #         ylim_max = axs_pf[ri][ci].get_ylim()[1]
    #         E_ref = 0.16
    #         axs_pf[ri][ci].plot(2*[E_ref],[0,1.1*ylim_max], '--k')
    #         axs_pf[ri][ci].text(x=E_ref+0.001, y=ylim_max/2, s=f'{E_ref} eV')
    #         axs_pf[ri][ci].set_ylim([0,ylim_max])

    # add heaviside demo
    first_plt_dct = comparison_lst[0]

    axs_pf[0][0].legend()
    axs_pf[1][0].legend()

    for spec_ax in [axs_pf[0][0], axs_pf[1][0]]:
        spec_ax.set_ylabel('Spectral Photon Flux Density, F$_\mathrm{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
        spec_ax.set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')

        # secax = spec_ax.secondary_xaxis('top', functions=(Eph_to_v, v_to_Ephs))
        # secax.set_xlabel('Wavenumber [cm$^{-1}$]')
        # secax.tick_params(axis='x', length=7)

        if log_atmdat:
            spec_ax.set_yscale('log')

    for Ndot_ax in [axs_pf[0][1], axs_pf[1][1]]:
        Ndot_ax.set_ylabel('Photon Flux Density, $\mathrm{\dot{N} \;[s^{-1}.m^{-2}]}$')
        Ndot_ax.set_xlabel('Bandgap, E$_\mathrm{g}$ [eV]')
        if log_atmdat:
            Ndot_ax.set_yscale('log')

        # secax = Ndot_ax.secondary_xaxis('top', functions=(Eph_to_v, v_to_Ephs))
        # secax.set_xlabel('Wavenumber, $\\tilde{v}$ [cm$^{-1}$]')

    for pf_ax in axs_pf.flatten():
        pf_ax.minorticks_on()
        pf_ax.yaxis.set_tick_params(which='minor', bottom=False)  # turn off minor ticks for y axis
        add_wn_ticks(pf_ax)
        og_ylim = pf_ax.get_ylim()
        pf_ax.set_ylim([0,og_ylim[1]])

    axs_pf[0][1].set_title('Photons emitted by TRD')
    axs_pf[1][1].set_title('Photons absorbed from environment')

if include_muopt_plots:
    if log_power:
        axs_Popt[0].set_yscale('log')
    axs_Popt[0].set_xlabel('Bandgap, E$_\mathrm{g}$ [eV]')
    axs_Popt[0].set_ylabel('Max Power Density [W.m$^{-2}$]')

    if use_cust_legend:
        axs_Popt[0].legend(handles=custom_muopt_legend, loc='upper right')

    axs_Popt[1].set_xlabel('Bandgap, E$_\mathrm{g}$ [eV]')
    axs_Popt[1].set_ylabel('V$_\mathrm{mpp}$ [V]')

    for ax in axs_Popt:
        if secondary_ticks == 'wavenumber':
            secax = ax.secondary_xaxis('top', functions=(Eph_to_v, v_to_Ephs))
            secax.set_xlabel('Wavenumber, $\\tilde{v}$ [cm$^{-1}$]')
            ax.minorticks_on()
            secax.minorticks_on()
        elif secondary_ticks == 'wavelength':
            add_wl_ticks(ax)


    if atmdat_background:
        axs_Popt[0].set_zorder(dwn_flux_yaxs.get_zorder()+1)
        axs_Popt[0].set_frame_on(False)
    # axs_Popt[0].set_xlim([0.0001,0.3])

if include_Eg_PD_scatter:
    axs_scatter.set_xlabel('Optimal Bandgap, E$_g$ [eV]')
    axs_scatter.set_ylabel('Max Power Density [W.m$^{-2}$]')
    axs_scatter.set_xlim([0.092, 0.102])
    axs_scatter.minorticks_on()
    axs_scatter.grid()

plt.show()
