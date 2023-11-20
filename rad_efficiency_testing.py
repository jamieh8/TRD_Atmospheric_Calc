import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.lines import Line2D



# Radiative efficiency testing:
comparison_lst = []
args_to_opt = ['mu']

# comparing different cwvs
# datsets = [{'cwv':10, 'Tc':296.724, 'colour':'deeppink'},
#            {'cwv':24, 'Tc':304.868, 'colour':'steelblue'},
#            {'cwv':54, 'Tc':303.512, 'colour':'seagreen'}]
#
# for ds in datsets:
#     cwv = ds['cwv']
#     atm_data = atmospheric_dataset(cwv=cwv)
#     Ephs = atm_data.photon_energies
#     line_format_dct = {'color': ds['colour'], 'linestyle': 'solid', 'alpha':0.75}
#     emitter_planck = planck_law_body(T=ds['Tc'], Ephs=Ephs)
#
#     comparison_lst += [{'label': f'cwv{cwv}', 'colour':ds['colour'],
#                         'cwv':cwv,
#                         'TRD in atm':TRD_in_atmosphere(emitter_planck, atm_data)}]

datasets = [
    {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o'},
    {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o'},
    {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o'},

    {'loc':'california', 'cwvstring':'low', 'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'symbol': 's'},
    {'loc':'california', 'cwvstring':'mid', 'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'symbol': 's'},
    {'loc':'california', 'cwvstring':'high', 'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'symbol': 's'},

    {'loc':'tamanrasset', 'cwvstring':'low', 'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'symbol':'^'},
    {'loc':'tamanrasset', 'cwvstring':'mid', 'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'symbol':'^'},
    {'loc':'tamanrasset', 'cwvstring':'high', 'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'symbol':'^'}
    ]

for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str)
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['color'], 'linestyle': 'solid'}
    scatter_format = {'c': ds['color'], 'marker': ds['symbol'], 'markersize':8}
    emitter_planck = planck_law_body(T=ds['Tskin'], Ephs=Ephs)
    comparison_lst += [{'label': f'{loc_str} {cwv_str}', 'colour':ds['color'], 'scatter format':scatter_format,
                        'cwv':ds['tcwv'],
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, atm_data)}]



emitter_planck = planck_law_body(T=300)
env_planck = planck_law_body(T=3)
comparison_lst += [{'label': f'3K', 'colour':'black','scatter format': {'c':'k', 'marker':'o', 'markersize':8},
                        'cwv':0,
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, env_planck)}]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)




fig_scatter, axs = plt.subplots(1,2,layout='tight', sharey='all', width_ratios=[1,2])
diodes = [{'Eg':0.094, 'eta':1, 'styleargs':{}}, {'Eg':0.094, 'eta':0.01, 'styleargs':{'fillstyle':'right'}}, {'Eg':0.25, 'eta':0.1, 'styleargs':{'mfc':'none', 'mew':2}}]
ax_scatter = axs[0]
for diode in diodes:
    for case in comparison_lst:
        opt_xs, opt_pd = get_best_pd(case['TRD in atm'],
                                     args_to_opt=args_to_opt,
                                     args_to_fix={'cutoff_angle':None, 'eta_ext':diode['eta'], 'Eg':diode['Eg'], 'consider_nonrad':True},
                                     alg=alg_powell)

        style_args = diode['styleargs']
        style_args.update(case['scatter format'])
        ax_scatter.plot(case['cwv'], opt_pd[0]*(-1), **style_args)


ax_scatter.set_yscale('log')
ax_scatter.set_xlabel('TCWV [mm]')
ax_scatter.set_ylabel('Power Density [W.m$^{-2}$]')
ax_scatter.minorticks_on()

reflines = [{'y':29, 'string':'yearly avg solar', 'xl':-5, 'xr':10, 'xt':10, 'arrow args':{'arrowstyle':'<-', 'color':'orangered'}, 'text args':{'color':'orangered', 'va':'center', 'ha':'left'}},
            {'y':51, 'string':'300K to 3K limit', 'xl':-5, 'xr':10, 'xt':10, 'arrow args':{'arrowstyle':'<-', 'color':'black'}, 'text args':{'color':'black', 'va':'center', 'ha':'left'}}]
for rl in reflines:
    ax_scatter.annotate(text='', xy=(rl['xr'],rl['y']), xytext=(rl['xl'],rl['y']), arrowprops=rl['arrow args'])
    ax_scatter.plot([0,80],2*[rl['y']], color='white', lw=1)
    ax_scatter.text(s=rl['string'], x=rl['xt'], y=rl['y'], **rl['text args'])



# Plot power magnitude examples
sample_area = 10  # [m-2]
power_magnitudes_guides = [{'label':'Average single person household \nin Australia for 1 day', 'PD':8*1e3 / (sample_area*24)},
                           {'label':'800W microwave for 15 mins', 'PD':800*0.25 / (sample_area*24)},
                           {'label':'10W LED bulb for 5h', 'PD':10*5 / (sample_area*24)},
                           {'label':'5W phone charger for 2h', 'PD':10 / (sample_area*24)},
                           {'label':'60W fridge for 24h', 'PD':60*24 / (sample_area*24)},]
ax_powerguide = axs[1]
for pguide in power_magnitudes_guides:
    ax_powerguide.annotate(pguide['label'], xy=(0,pguide['PD']), xytext=(0.1, pguide['PD']), arrowprops=dict(arrowstyle="->"), va='center')
    # ax_powerguide.text(s=pguide['label'], x=0.05, y=pguide['PD'], ha='left', va='bottom')
    print(pguide['label'] + ' ' + str(pguide['PD']*24))
ax_powerguide.set_title(f'Over 24h, {sample_area} m$^2$ can power:')
ax_powerguide.axis('off')

# Add legend in second plot
custom_muopt_legend = []
for diode in diodes:
    Eg, eta = diode['Eg'], diode['eta']
    custom_muopt_legend += [Line2D([0],[0], linestyle='none', **diode['styleargs'],
                                   label='E$_\mathrm{g}$=' + f'{Eg} eV' + ', $\eta_\mathrm{ext}=$' + f'{eta*100:.0f}%')]

ax_powerguide.legend(handles=custom_muopt_legend, loc='lower left')

# ax_scatter.grid()
ax_scatter.set_xlim([-5,75])
ax_scatter.set_ylim([1e-4, 1e2])
# ax_powerguide.set_ylim([1e-4, 1e2])
# sec_yaxs = ax_scatter.secondary_yaxis('right', functions=(lambda x: x*24, lambda x: x/24))
# sec_yaxs.set_ylabel('Energy Density over 24h [Wh.m$^{-2}$]')
# axs_scatter.set_title('$\eta_\mathrm{ext}$=10$^{-2}$, $E_\mathrm{g}=0.1$ eV')


# ------------------------ Sweep eta_ext ------------------------
# Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
# cutoff_angle = None
#
# cwv, Tc = 10, 296.724
# emitter = planck_law_body(Tc, Ephs)
# Eg = 0.1
# schemes_to_plot = []
#
# # for Te, linestyle in zip([3, 210], ['-','--']):
# #     atm_bb = planck_law_body(Te, Ephs)
# #     combined_atmBB = TRD_in_atmosphere(emitter, atm_bb)
# #     schemes_to_plot += [{'label': f'{Te:.0f} env', 'colour': 'black', 'linestyle': linestyle,
# #                  'Eg': Eg, 'TRD in atm': combined_atmBB}]
#
# for cwv, cwv_col in zip([10],['deeppink']): #zip([10,24,54],['deeppink', 'steelblue', 'seagreen']):
#     atm_dataset = atmospheric_dataset(cwv)
#     combined_cwv = TRD_in_atmosphere(emitter, atm_dataset)
#     schemes_to_plot += [{'label':f'cwv{cwv}, '+'$\eta_{rad}$=', 'colour':cwv_col, 'linestyle':'-',
#                  'Eg':Eg, 'TRD in atm': combined_cwv}]

# fig_rs, axs_rs = plt.subplots(1, 3, layout='tight')
# log_eta = True
# if log_eta:
#     rad_eff_sweep = np.logspace(-4, 0, num=10, base=10)
# else:
#     rad_eff_sweep = np.linspace(0.01,1,10)
#
# alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
# args_to_opt = ['mu', 'Eg']
#
#
# relative_change = True
#
# for scheme in schemes_to_plot:
#     TRD_in_atm_obj = scheme['TRD in atm']
#     i_colour = scheme['colour']
#     # Sweep rad efficiency
#     maxPds = []
#     opt_vals = []
#     for eta in rad_eff_sweep:
#         opt_xs, opt_pd = get_best_pd(TRD_in_atm_obj, args_to_opt=args_to_opt,
#                                      args_to_fix={'cutoff_angle':90, 'eta_ext':eta, 'consider_nonrad':True}, alg=alg_powell)
#         opt_vals += [list(opt_xs.values())]
#         maxPds += [opt_pd[0]]
#
#     if relative_change:
#         ys1 = maxPds/maxPds[-1]
#         # ys2 = Vmpps/Vmpps[-1]
#     else:
#         ys1 = maxPds
#
#     if log_eta:
#         axs_rs[0].loglog(rad_eff_sweep, ys1, scheme['linestyle'], label = scheme['label'], color=i_colour)  # current
#     else:
#         axs_rs[0].plot(rad_eff_sweep, ys1, scheme['linestyle'], label=scheme['label'], color=i_colour)  # current
#
#     # Secondary plot showing x values corresponding to best P.D.
#     opt_vals_t = np.transpose(np.array(opt_vals))
#     for ir, row in enumerate(opt_vals_t):
#         if log_eta:
#             axs_rs[ir+1].loglog(rad_eff_sweep, np.abs(row), scheme['linestyle'], c=scheme['colour'], label=scheme['label'])
#         else:
#             axs_rs[ir + 1].loglog(rad_eff_sweep, row, scheme['linestyle'], c=scheme['colour'],label=scheme['label'])
#
# for ax in axs_rs:
#     ax.set_xlabel('$\eta_{ext}$')
#
# translate_to_label = {'Eg':'Bandgap, E$_\mathrm{g}$ [eV]', 'mu':'V$_\mathrm{mpp}$ [V]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}
# if 'mu' in args_to_opt:
#     args_to_opt.remove('mu')
#     args_to_opt += ['mu']
# for opi, opt_label in enumerate(args_to_opt):
#     axs_rs[opi+1].set_ylabel(translate_to_label[opt_label])
#
# if relative_change:
#     axs_rs[0].set_ylabel('P / P($\eta$=1)')
# else:
#     axs_rs[0].set_ylabel('Max Power Density [W.m$^{-2}$]')


plt.show()
