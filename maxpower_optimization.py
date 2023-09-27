import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar
from matplotlib.colors import Normalize, LogNorm



# ----------------- Line/Sweep Optimization Plotter -----------------

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

fix_Eg = False

# Sweep cutoff angle, optimize over Eg and mu
# args_to_opt = ['Eg', 'mu']
# x_vals_sweep = np.arange(10,95,5)
# arg_sweep_str = 'cutoff_angle'
# arg_fix_extra = {}
# xlabel_str = 'Cutoff Angle [$\circ$]'

# Sweep Eg, fix cutoff angle, optimize over mu
# args_to_opt = ['mu']
# x_vals_sweep = np.linspace(0.062,0.2,20)  # bandgaps, in [eV]
# arg_sweep_str = 'Eg'
# arg_fix_extra = {'cutoff_angle':None, 'eta_ext':0.1, 'consider_nonrad':True}
# xlabel_str = 'Bandgap, E$_g$ [eV]'

# Sweep Eg, optimize over mu and cutoff angle
args_to_opt = ['mu', 'cutoff_angle']
x_vals_sweep = np.linspace(0.062,0.2,20)  # bandgaps, in [eV]
arg_sweep_str = 'Eg'
arg_fix_extra = {'eta_ext':1, 'consider_nonrad':False}
xlabel_str = 'Bandgap, E$_g$ [eV]'

# Sweep cutoff angle, optimize over mu (with Eg defined per case)
# args_to_opt = ['mu']
# x_vals_sweep = np.arange(10,95,5)
# arg_sweep_str = 'cutoff_angle'
# xlabel_str = 'Cutoff Angle, $\\theta_c$ [$\circ$]'
# fix_Eg = True

to_plot = []


# Setup cases to plot
for AD_dict in [
    {'cwv':10, 'Tc':296.724, 'col':'deeppink', 'Eg':0.094},
    {'cwv':24, 'Tc':304.868, 'col':'steelblue', 'Eg':0.094},
    {'cwv':54, 'Tc':303.512, 'col':'seagreen', 'Eg':0.1}
]:
    Eg = AD_dict['Eg']
    cwv = AD_dict['cwv']
    atm_data = atmospheric_dataset(cwv)
    TRD = planck_law_body(T=AD_dict['Tc'], Ephs=atm_data.photon_energies)
    comb_TRDenv = TRD_in_atmosphere(TRD, atm_data)
    if fix_Eg:
        arg_fix_extra.update({'Eg':Eg})  #<--- comment out if Eg is to be optimized
    to_plot += [{'TRD in env':comb_TRDenv, 'label':f'cwv{cwv}', 'colour':AD_dict['col'], 'arg_fix_extra': arg_fix_extra}]

# TRD_300K = planck_law_body(T=300, Ephs=Ephs)
# for Te, colT in zip([3,200],['black','dimgrey']):
#     env_T = planck_law_body(T=Te, Ephs=Ephs)
#     comb_TRDenv = TRD_in_atmosphere(TRD_300K, env_T)
#     to_plot += [{'TRD in env':comb_TRDenv, 'label':f'300K/{Te}K, E$_g$=0.02', 'colour':colT, 'arg_fix_extra': {'Eg':0.02}}]


# Define subplots based on optimization
fig, axs = plt.subplots(1,len(args_to_opt)+1, layout='tight')
relative_change = False

for case in to_plot:
    print(case['label'])
    combined_trd_env = case['TRD in env']
    arg_fix_extra = case['arg_fix_extra']

    max_pds = []
    opt_vals = []
    max_pds_90, mus_90 = [], []
    for xval in x_vals_sweep:
        arg_f = {arg_sweep_str: xval}
        if arg_fix_extra != None:
            arg_f.update(arg_fix_extra)
        opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg = alg_de)
        max_pds += [opt_pd[0]]
        opt_vals += [[opt_xs.get(str_lbl) for str_lbl in args_to_opt]]  # ensures order of values in opt_vals matches order of args_to_opt.!

        print(f'optimized xs: {opt_xs}')
        print(f'fixed: {arg_f}')
        print(f'max PD: {opt_pd[0]}')

        optres_90, pd_90 = get_best_pd(combined_trd_env, args_to_opt=['mu'], args_to_fix={arg_sweep_str: xval, 'cutoff_angle':None}, alg = alg_powell)
        max_pds_90 += [pd_90[0]]
        mus_90 += [optres_90['mu']]

    if relative_change:
        pd_ys = max_pds/max_pds[-1]
    else:
        pd_ys = max_pds
    axs[0].plot(x_vals_sweep, pd_ys, c=case['colour'], label=case['label'])

    # Secondary plot showing x values corresponding to best P.D.
    opt_vals_t = np.transpose(np.array(opt_vals))
    for ir, row in enumerate(opt_vals_t):
        axs[ir+1].plot(x_vals_sweep, row, c=case['colour'], label=case['label'])

    # reference plots for cutoff 90
    axs[0].plot(x_vals_sweep, max_pds_90, '--', c=case['colour'])
    axs[args_to_opt.index('cutoff_angle')+1].plot(x_vals_sweep, len(x_vals_sweep)*[90],'--', c=case['colour'])
    axs[args_to_opt.index('mu')+1].plot(x_vals_sweep, mus_90, '--', c=case['colour'])

pd_ax = axs[0]
if relative_change:
    pd_ax.set_ylabel(f'P / P({x_vals_sweep[-1]:.0f})')
else:
    pd_ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
pd_ax.set_xlabel(xlabel_str)
pd_ax.legend()

translate_to_label = {'Eg':'Bandgap, E$_g$ [eV]', 'mu':'V$_{mpp}$ [V]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}
for opi, opt_label in enumerate(args_to_opt):
    axs[opi+1].set_ylabel(translate_to_label[opt_label])
    axs[opi+1].set_xlabel(xlabel_str)



# ----------------- Heatmap ----------------- #
# fig, axs = plt.subplots(1,1)
# # Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
#
# # mus = np.linspace(-0.1, 0, 100)
# # cutoff_angle = 70
#
# # Cutoff Angle & Bandgap. Opt over mu
# x_sweep = np.linspace(10,90,21)
# x_id_str = 'cutoff_angle'
# Eg_start, Eg_end, Eg_count = 0.062, 0.15, 20
# y_sweep = np.linspace(Eg_start, Eg_end, Eg_count)
# y_id_str = 'Eg'
# args_to_opt= ['mu']
# arg_fix_extra = {'eta_ext':1, 'consider_nonrad':False}
# commercial_diode_ref = False
# norm_str = 'default'
#
# # Rad efficiency & Bandgap. Opt over mu.
# # x_sweep = np.linspace(0.01,1, 40)
# # x_sweep = np.logspace(-4,0,num=40,base=10)
# # x_id_str = 'eta_ext'
# # x_log = True
# # Eg_start = 0.062
# # y_sweep = np.linspace(Eg_start, 0.2, 40)
# # y_id_str = 'Eg'
# # args_to_opt = ['mu']
# # arg_fix_extra = {'cutoff_angle':None, 'consider_nonrad':True}
# # commercial_diode_ref = True
# # norm_str = 'log power'
#
# alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# # alg = pg.de(gen=50, ftol=1e-5)
#
#
# # cwv, Tc = 10, 296.724
# # cwv, Tc = 24, 304.868
# cwv, Tc = 54, 303.512
# atm_data = atmospheric_dataset(cwv)
# case_label = f'cwv{cwv}'
# # filename = f'PD_cwv{cwv}_etaextlog_-4_0_40_Egs_{Eg_start}_02_40.csv'
# filename = f'PD_cwv{cwv}_cutoffangle_10_90_21_Egs_{Eg_start}_{Eg_end}_{Eg_count}.csv'
#
# # Tc, Te = 300, 3
# # atm_data = planck_law_body(Te, Ephs)
# # case_label = f'T$_c$ {Tc}K, T$_e$ {Te}K'
# # filename = f'PD_Te3K_Tc300K_etaextlog_-4_0_40_Egs_{Eg_start}_02_40.csv'
#
# emitter_planck = planck_law_body(Tc, atm_data.photon_energies)
# combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data)
#
# relative_to_etaext1 = False
#
# # Generate new data
# # pds_2d = []
# # for y in y_sweep:
# #     row = []
# #     for x in x_sweep:
# #     # for mu in mus:
# #     #     if mu < -Eg:
# #     #         pd = np.nan
# #     #     else:
# #     #         pd = combined_trd_env.power_density(mu, Eg, cutoff_angle=cutoff_angle)
# #         arg_f = {x_id_str:x, y_id_str:y}
# #         if arg_fix_extra != None:
# #             arg_f.update(arg_fix_extra)
# #         # best_pd = [1]
# #         # while best_pd[0] > 0:
# #         best_xs, best_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg=alg)
# #         # print(best_pd[0])
# #         row += [best_pd[0]]
# #     pds_2d += [row]
# #
# # pds_2d = np.asarray(pds_2d)
# #
# # # Save file
# # np.savetxt(filename, pds_2d, delimiter = ',')
#
# # Import existing data
# pds_2d = np.loadtxt(filename, delimiter=',', dtype=float)
#
#
# # Test optimizer
# alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# alg = pg.de(gen=50, ftol=1e-4)
# # Egs_opt = []
# # Vs_opt = []
# # pds_opt = []
# for i in range(5):
#     opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=[y_id_str]+args_to_opt+[x_id_str], args_to_fix=arg_fix_extra, alg = alg)
#     plt.plot(opt_xs[x_id_str], opt_xs[y_id_str], 'x', c='black')
#     print(opt_pd[0])
#     print(opt_xs)
#     # Egs_opt += [opt_res[0][1]]
#     # Vs_opt += [opt_res[0][0]]
#     # pds_opt += [opt_res[1][0]]
#
# # for str, array in zip(['Egs','Vs','PDs'], [Egs_opt, Vs_opt, pds_opt]):
# #     print(str)
# #     for elem in array:
# #         print(elem)
#
#
# # Plot heatmap
# h_ax = axs
#
# if norm_str == 'mid0':
#     # used for Eg vs V, to identify areas of TRD operation
#     norm_0mid = Normalize(vmin=-10, vmax=10)
#     hmap = h_ax.pcolor(x_sweep, y_sweep, pds_2d, cmap='bwr', norm=norm_0mid, shading='nearest')  # heatmap
#     add_loglvls = False
#
# elif norm_str == 'log power':
#     # used for Eg vs rad efficiency, plotted logarithmically
#     if relative_to_etaext1:
#         # Relative to efficiency = 1 (last x)
#         pds_etaext1 = pds_2d[:, -1]
#         pds_relative = np.divide(pds_2d, pds_etaext1[:, np.newaxis])  # newaxis required to divide along correct axis
#         C_2darray = pds_relative
#         cb_label = r'P$_\mathrm{max}$ / P$_\mathrm{max}$ ($\eta_\mathrm{ext}$ = 1)'
#
#     else:
#         # Not relative - convert to + values so log scale works
#         C_2darray = (-1)*np.array(pds_2d)
#         cb_label = r'Power Density [W.m$^{-2}$]'
#     min_a = np.min(C_2darray)
#     max_a = np.max(C_2darray)
#     norm_log = LogNorm(vmin=min_a, vmax=max_a)
#     hmap = h_ax.pcolormesh(x_sweep, y_sweep, C_2darray, cmap='inferno', shading='gouraud', norm=norm_log)
#     add_loglvls = True
#
# else:
#     # default.. linear scale, no normalization.
#     hmap = h_ax.pcolormesh(x_sweep, y_sweep, (-1)*pds_2d, cmap='inferno', shading='gouraud')
#     add_loglvls = False
#     cb_label = r'Power Density [W.m$^{-2}$]'
#
#
# # Colorbar formatting
# cbar = plt.colorbar(hmap)
# cbar.ax.tick_params(labelsize=10)
# # cbar.ax.set_ylabel(r'Spectral Photon Flux [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]')
# cbar.ax.set_ylabel(cb_label)
#
# if add_loglvls:
#     log_lvls = 10. ** np.arange(start=-7, stop=2)
#     cntr = h_ax.contour(x_sweep, y_sweep, C_2darray, levels=log_lvls, colors='black')
#     for lvl in log_lvls:
#         cbar.ax.plot([0, 1], 2 * [lvl], '-k')
#
#
# # adding commercial diodes for ref
# if commercial_diode_ref:
#     diodes = [{'eta_ext':10**(-2), 'Eg':0.1}, {'eta_ext':10**(-2), 'Eg':0.175}, {'eta_ext':10**(-2), 'Eg':0.2}]
#     for diode in diodes:
#         plt.plot(diode['eta_ext'], diode['Eg'], 'or')
#         ind_x = x_sweep.searchsorted(diode['eta_ext'])
#         ind_y = y_sweep.searchsorted(diode['Eg'])
#         cbar.ax.plot([0.5], [C_2darray[ind_y][ind_x]], 'or')
#
# # Labeling axes
# translate_to_label = {'Eg':'Bandgap, E$_\mathrm{g}$ [eV]', 'mu':'V$_\mathrm{mpp}$ [V]',
#                       'cutoff_angle': 'Cutoff Angle, $\\theta_\mathrm{c}$ [$\circ$]', 'eta_ext':'Rad Efficiency, $\eta_\mathrm{ext}$'}
# h_ax.set_xlabel(translate_to_label[x_id_str])
# h_ax.set_ylabel(translate_to_label[y_id_str])
# h_ax.set_xlim([x_sweep[0], x_sweep[-1]])
#
#
# plt.title(f'Optimized over $\mu$, {case_label}')

plt.show()