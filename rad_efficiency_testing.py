import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *




# Radiative efficiency testing:
comparison_lst = []
args_to_opt = ['mu']
# comparing different cwvs
Egs_AD = np.arange(0.062, 0.3, 0.002)
datsets = [{'cwv':10, 'Tc':296.724, 'colour':'deeppink'},
           {'cwv':24, 'Tc':304.868, 'colour':'steelblue'},
           {'cwv':54, 'Tc':303.512, 'colour':'seagreen'}]
for ds in datsets:
    cwv = ds['cwv']
    atm_data = atmospheric_dataset(cwv=cwv)
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['colour'], 'linestyle': 'solid', 'alpha':0.75}
    emitter_planck = planck_law_body(T=ds['Tc'], Ephs=Ephs)

    comparison_lst += [{'label': f'cwv{cwv}', 'colour':ds['colour'],
                        'cwv':cwv,
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, atm_data)}]

emitter_planck = planck_law_body(T=300)
env_planck = planck_law_body(T=3)
comparison_lst += [{'label': f'3K', 'colour':'black',
                        'cwv':0,
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, env_planck)}]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)




fig_scatter, axs_scatter = plt.subplots(1,1,layout='tight')
Eg = 0.1
for eta, symbol in zip([1e-1, 1e-2], ['o','s']):
    for case in comparison_lst:
        opt_xs, opt_pd = get_best_pd(case['TRD in atm'],
                                     args_to_opt=args_to_opt,
                                     args_to_fix={'cutoff_angle':None, 'eta_ext':eta, 'Eg':Eg, 'consider_nonrad':True},
                                     alg=alg_powell)
        axs_scatter.plot(case['cwv'], opt_pd[0]*(-1), symbol, c=case['colour'])

axs_scatter.set_yscale('log')
axs_scatter.set_xlabel('cwv')
axs_scatter.set_ylabel('Power Density [W.m$^{-2}$]')
axs_scatter.minorticks_on()

reflines = [{'y':1e6/24, 'string':'daily average solar'}, {'y':51, 'string':'300K to 3K limit'}]
for rl in reflines:
    axs_scatter.plot([0,60],2*[rl['y']], '-k', lw=1)
    axs_scatter.text(s=rl['string'],x=25, y=rl['y'], ha='center', va='bottom')


axs_scatter.set_xlim([0,55])
sec_yaxs = axs_scatter.secondary_yaxis('right', functions=(lambda x: x*24/1e3, lambda x: x/24 * 1e3))

sec_yaxs.set_ylabel('Energy Density over 24h [Wh.m$^{-2}$]')
axs_scatter.set_title('$\eta_\mathrm{ext}$=10$^{-2}$, $E_\mathrm{g}=0.1$ eV')

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
