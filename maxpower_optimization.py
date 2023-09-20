import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar
from matplotlib.colors import Normalize, LogNorm



# ----------------- Line/Sweep Optimization Plotter -----------------

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

# Sweep cutoff angle, optimize over Eg and mu
# args_to_opt = ['Eg', 'mu']
# x_vals_sweep = np.arange(10,95,10)
# arg_sweep_str = 'cutoff_angle'
# arg_fix_extra = {}
# xlabel_str = 'Cutoff Angle [$\circ$]'

# Sweep Eg, fix cutoff angle, optimize over mu
args_to_opt = ['mu']
x_vals_sweep = np.linspace(0.062,0.2,20)  # bandgaps, in [eV]
arg_sweep_str = 'Eg'
arg_fix_extra = {'cutoff_angle':None, 'eta_ext':0.1, 'consider_nonrad':True}
xlabel_str = 'Bandgap, E$_g$ [eV]'

# Sweep Eg, optimize over mu and cutoff angle
# args_to_opt = ['mu', 'cutoff_angle']
# x_vals_sweep = np.linspace(0.062,0.15,20)  # bandgaps, in [eV]
# arg_sweep_str = 'Eg'
# arg_fix_extra = {}
# xlabel_str = 'Bandgap, E$_g$ [eV]'

# Sweep cutoff angle, optimize over mu (with Eg defined per case)
# args_to_opt = ['mu']
# x_vals_sweep = np.arange(10,95,5)
# arg_sweep_str = 'cutoff_angle'
# xlabel_str = 'Cutoff Angle, $\\theta_c$ [$\circ$]'

to_plot = []


# Setup cases to plot
for AD_dict in [
    {'cwv':10, 'Tc':296.724, 'col':'deeppink', 'Eg':0.094}#,
    # {'cwv':24, 'Tc':304.868, 'col':'steelblue', 'Eg':0.094},
    # {'cwv':54, 'Tc':303.512, 'col':'seagreen', 'Eg':0.1}
]:
    Eg = AD_dict['Eg']
    cwv = AD_dict['cwv']
    TRD = planck_law_body(T=AD_dict['Tc'], Ephs=Ephs)
    atm_data = atmospheric_dataset(cwv)
    comb_TRDenv = TRD_in_atmosphere(TRD, atm_data)
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
    combined_trd_env = case['TRD in env']
    arg_fix_extra = case['arg_fix_extra']

    max_pds = []
    opt_vals = []
    max_pds_90, mus_90 = [], []
    for xval in x_vals_sweep:
        arg_f = {arg_sweep_str: xval}
        if arg_fix_extra != None:
            arg_f.update(arg_fix_extra)
        opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg = alg_powell)
        max_pds += [opt_pd[0]]
        opt_vals += [list(opt_xs.values())]

        # mu_90, pd_90 = get_best_pd(combined_trd_env, args_to_opt=['mu'], args_to_fix={arg_sweep_str: xval, 'cutoff_angle':90}, alg = alg_powell)
        # max_pds_90 += [pd_90]
        # mus_90 += [mu_90]

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
    # axs[0].plot(x_vals_sweep, max_pds_90, '--', c=case['colour'])
    # axs[1].plot(x_vals_sweep, mus_90, '--', c=case['colour'])
    # axs[2].plot(x_vals_sweep, len(x_vals_sweep)*[90],'--', c=case['colour'])

pd_ax = axs[0]
if relative_change:
    pd_ax.set_ylabel(f'P / P({x_vals_sweep[-1]:.0f})')
else:
    pd_ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
pd_ax.set_xlabel(xlabel_str)
pd_ax.legend()

translate_to_label = {'Eg':'Bandgap, E$_g$ [eV]', 'mu':'V$_{mpp}$ [V]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}
if 'mu' in args_to_opt:
    args_to_opt.remove('mu')
    args_to_opt += ['mu']
for opi, opt_label in enumerate(args_to_opt):
    axs[opi+1].set_ylabel(translate_to_label[opt_label])
    axs[opi+1].set_xlabel(xlabel_str)



# ----------------- Heatmap ----------------- #
# fig, axs = plt.subplots(1,1)
# Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
# # mus = np.linspace(-0.1, 0, 100)
# # cutoff_angle = 70
#
# # Cutoff Angle & Bandgap. Opt over mu
# # x_sweep = np.linspace(10,90,20)
# # x_id_str = 'cutoff_angle'
# # y_sweep = np.linspace(0.001, 0.15, 40)
# # y_id_str = 'Eg'
# # args_to_opt= ['mu']
#
# # rad efficiency & bandgap. opt over mu.
# # x_sweep = np.linspace(0.01,1, 40)
# x_sweep = np.logspace(-4,0,num=40,base=10)
# x_id_str = 'eta_ext'
# x_log = True
# Eg_start = 0.062
# y_sweep = np.linspace(Eg_start, 0.2, 40)
# y_id_str = 'Eg'
# args_to_opt = ['mu']
# arg_fix_extra = {'cutoff_angle':None, 'consider_nonrad':True}
#
# # alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# alg = pg.de(gen=50, ftol=1e-5)
#
#
# # cwv, Tc = 10, 296.724
# # cwv, Tc = 24, 304.868
# cwv, Tc = 54, 303.512
# atm_data = atmospheric_dataset(cwv)
# case_label = f'cwv{cwv}'
# filename = f'PD_cwv{cwv}_etaextlog_-4_0_40_Egs_{Eg_start}_02_40.csv'
#
# # Tc, Te = 300, 3
# # atm_data = planck_law_body(Te, Ephs)
# # case_label = f'T$_c$ {Tc}K, T$_e$ {Te}K'
# # filename = f'PD_Te3K_Tc300K_etaextlog_-4_0_40_Egs_{Eg_start}_02_40.csv'
#
# emitter_planck = planck_law_body(Tc, Ephs)
# combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data)
#
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
# #         best_pd = [1]
# #         while best_pd[0] > 0:
# #             best_xs, best_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg=alg)
# #         row += [best_pd[0]]
# #     pds_2d += [row]
#
#
# # Save file
# # np.savetxt(filename, np.asarray(pds_2d), delimiter = ',')
#
# # Import existing data
# pds_2d = np.loadtxt(filename, delimiter=',', dtype=float)
#
#
# # Test optimizer
# # alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# # # alg = pg.de(gen=50, ftol=1e-5)
# # Egs_opt = []
# # Vs_opt = []
# # pds_opt = []
# # for i in range(5):
# #     opt_res = get_best_pd(combined_trd_env, args_to_opt={'Eg':0, 'mu':0}, args_to_fix={'cutoff_angle':cutoff_angle}, alg = alg)
# #     plt.plot([opt_res[0][1]], [opt_res[0][0]], 'x', c='black')
# #
# #     Egs_opt += [opt_res[0][1]]
# #     Vs_opt += [opt_res[0][0]]
# #     pds_opt += [opt_res[1][0]]
# #
# # for str, array in zip(['Egs','Vs','PDs'], [Egs_opt, Vs_opt, pds_opt]):
# #     print(str)
# #     for elem in array:
# #         print(elem)
#
# # Plot heatmap
# h_ax = axs
#
#
# # norm_0mid = Normalize(vmin=-10, vmax=10)
# # hmap = h_ax.pcolor(mus, Egs, pds_2d, cmap='bwr', norm=norm_0mid, shading='nearest')  # heatmap
# # h_ax.set_xlabel('V [V] / $\mu$ [eV]')
#
# pds_2d = (-1)*np.array(pds_2d)
# min_pd = np.min(pds_2d)
# max_pd = np.max(pds_2d)
# print(f'{min_pd}, {max_pd}')
# norm_log = LogNorm(vmin=min_pd, vmax=max_pd)
# hmap = h_ax.pcolormesh(x_sweep, y_sweep, pds_2d, cmap='inferno', shading='gouraud', norm=norm_log)
# log_lvls = 10. ** np.arange(start=-7, stop=2)
# cntr = h_ax.contour(x_sweep, y_sweep, pds_2d, levels=log_lvls, colors='black')
#
# translate_to_label = {'Eg':'Bandgap, E$_\mathrm{g}$ [eV]', 'mu':'V$_\mathrm{mpp}$ [V]',
#                       'cutoff_angle': 'Cutoff Angle, $\\theta$ [$\circ$]', 'eta_ext':'Rad Efficiency, $\eta_\mathrm{ext}$'}
# h_ax.set_xlabel(translate_to_label[x_id_str])
# h_ax.set_ylabel(translate_to_label[y_id_str])
# h_ax.set_xscale('log')
# h_ax.set_xlim([x_sweep[0], x_sweep[-1]])
#
# cbar = plt.colorbar(hmap)
# for lvl in log_lvls:
#     cbar.ax.plot([0,1], 2*[lvl], '-k')
# cbar.ax.tick_params(labelsize=10)
# # cbar.ax.set_ylabel(r'Spectral Photon Flux [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]')
# cbar.ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
#
# plt.title(f'Optimized over $\mu$, {case_label}')

plt.show()