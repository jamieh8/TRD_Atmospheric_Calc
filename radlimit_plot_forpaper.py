import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def Eph_to_v(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def v_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')


comparison_lst = []

# comparing new datasets:
angle_array = [0]
Egs_AD = np.arange(0.0125, 0.3, 0.002)
datasets_telfer = get_dataset_list()[0:3]

for ds in datasets_telfer:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str, Tskin=ds['Tskin'], date='23dec')
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['color'], 'linestyle': 'solid'}
    scatter_format = {'c': ds['color'], 'marker': ds['symbol'], 'markersize':8}
    emitter_planck = planck_law_body(T=ds['Tskin'], Ephs=Ephs)
    comparison_lst += [{'label': f'{loc_str.capitalize()} {cwv_str}', 'line format':line_format_dct, 'scatter format':scatter_format,
                        'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
                         'cutoff angle':None, 'use diffusivity approx':True, 'Egs':Egs_AD}]


# comparing blackbody environments
Ephs = np.arange(1e-6, 0.31, 0.0001)
Egs_bb = np.arange(0.0001,0.3,0.002)
emitter_planck_300 = planck_law_body(T=300, Ephs=Ephs)
Tsets = [{'Tc':3, 'colour':'black'}]
#         #, {'Tc':200, 'colour':'navy'}, {'Tc':270, 'colour':'blueviolet'}, {'Tc':290, 'colour':'mediumorchid'}]
for dataset_entry in comparison_lst:
    Teffective = dataset_entry['atmospheric dataset'].effective_skytemp(300)
    Tsets += [{'Tc':Teffective, 'colour':dataset_entry['line format']['color']}]

for Ts in Tsets:
    Tc = Ts['Tc']
    bb_env = planck_law_body(Tc, Ephs)
    line_format_dct = {'color': Ts['colour'], 'linestyle': 'dashed'}
    comparison_lst += [{'label': 'T$_\mathrm{atm}$ = '+f'{Tc:.5g}K', 'line format':line_format_dct,
                        'atmospheric dataset': bb_env, 'emitter body': emitter_planck_300,
                        'cutoff angle': None, 'use diffusivity approx': True, 'Egs':Egs_bb}]

comparison_lst[-1].update({'label position':'below'})  # required for Tatm Telfer high, to accomodate Telfer mid


custom_muopt_legend = [Line2D([0],[0], color = 'k', linestyle='solid', label='LBLTRM modelling'),
                       Line2D([0],[0], color = 'k', linestyle='dashed', label='Effective temperature'),
                       Patch(facecolor='darkorange', label='Telfer low'),
                       Patch(facecolor='darkviolet', label='Telfer mid'),
                       Patch(facecolor='mediumseagreen', label='Telfer high'),
                       Patch(facecolor='black', label='3K BB')]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

fig, axs = plt.subplots(1, 2, width_ratios=[3,1])

ref_dataset = atmospheric_dataset_new('low', 'telfer', Tskin=301.56, date='24oct')
downwelling_photflux = ref_dataset.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV',
                                                               col_name='downwelling_flux')
dwn_flux_yaxs = axs[0].twinx()
dwn_flux_yaxs.plot(ref_dataset.photon_energies, downwelling_photflux, c='lightgrey')
dwn_flux_yaxs.set_ylabel('Spectral Photon Flux Density, $\mathrm{F_{ph} \; [s^{-1}.m^{-2}/eV]}$', color='lightgrey')
dwn_flux_yaxs.tick_params(axis='y', labelcolor='lightgrey')
# dwn_flux_yaxs.text(s='Telfer low', x=0.09, y=0.5 * 1e23, ha='right', color='lightgrey')
dwn_flux_yaxs.set_yscale('log')

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

    combined_obj = TRD_in_atmosphere(emitter, atm_data)
    # TRD performance metrics
    maxPs = []
    Vmpps = []
    for Eg in Egs:
        mu_opt = combined_obj.optimize_mu(Eg, cutoff_angle)
        maxPs += [mu_opt['max power']]
        Vmpps += [mu_opt['Vmpp']]

    maxPs = (-1)*np.array(maxPs)
    axs[0].plot(Egs, maxPs, **style_args, label=sample_dct['label'])

    opt_xs, opt_pd = get_best_pd(combined_obj, args_to_opt=['Eg','mu'],
                                args_to_fix={'cutoff_angle':None, 'consider_nonrad':False, 'eta_ext':1}, alg=alg_powell)

    pd = opt_pd[0] * (-1)
    y_offset = 0.15*pd
    axs[0].plot(opt_xs['Eg'], pd, 'o', **style_args)

axs[0].set_yscale('log')
axs[0].set_xlabel('Bandgap, E$_\mathrm{g}$ [eV]')
axs[0].set_ylabel('Max Power Density [W.m$^{-2}$]')


axs[0].legend(handles=custom_muopt_legend, loc='upper right')

add_wl_ticks(axs[0])
axs[0].set_xlim([0.009,0.31])
axs[0].set_zorder(dwn_flux_yaxs.get_zorder()+1)
axs[0].set_frame_on(False)



# For scatter plot
axs_scatter = axs[1]
scatter_lst = []
angle_array = [0]
Egs_AD = np.arange(0.0125, 0.3, 0.002)

datasets_all = get_dataset_list()

for ds in datasets_all:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str, Tskin=ds['Tskin'], date='23dec')
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['color'], 'linestyle': 'solid'}
    scatter_format = {'c': ds['color'], 'marker': ds['symbol'], 'markersize':8}
    emitter_planck = planck_law_body(T=ds['Tskin'], Ephs=Ephs)
    scatter_lst += [{'label': f'{loc_str.capitalize()} {cwv_str}', 'line format':line_format_dct, 'scatter format':scatter_format,
                        'atmospheric dataset':atm_data, 'emitter body':emitter_planck,
                         'cutoff angle':None, 'use diffusivity approx':True, 'Egs':Egs_AD}]

for sample_dct in scatter_lst:
    atm_data = sample_dct['atmospheric dataset']
    emitter = sample_dct['emitter body']
    combined_obj = TRD_in_atmosphere(emitter, atm_data)
    opt_xs, opt_pd = get_best_pd(combined_obj, args_to_opt=['Eg', 'mu'],
                                 args_to_fix={'cutoff_angle': None, 'consider_nonrad': False, 'eta_ext': 1}, alg=alg_de)
    pd = opt_pd[0] * (-1)

    axs_scatter.plot(opt_xs['Eg'], pd, **sample_dct['scatter format'])

axs_scatter.set_xlabel('Optimal Bandgap, E$_\mathrm{g}$ [eV]')
axs_scatter.set_ylabel('Max Power Density [W.m$^{-2}$]')
axs_scatter.set_xlim([0.092, 0.102])
axs_scatter.minorticks_on()
axs_scatter.grid()


axs[0].set_title('a)', loc='left', fontsize=15, fontweight='bold', y=1.05, x=-0.05)
axs_scatter.set_title('b)', loc='left', fontsize=15, fontweight='bold', y=1.05, x=-0.15)


plt.show()