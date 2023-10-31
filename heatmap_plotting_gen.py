import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker

# ----------------- Heatmap ----------------- #
fig, axs = plt.subplots(1,1)
# Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

# mus = np.linspace(-0.1, 0, 100)
# cutoff_angle = 70

# Cutoff Angle & Bandgap. Opt over mu
# x_sweep = np.linspace(10,90,21)
# x_id_str = 'cutoff_angle'
# Eg_start, Eg_end, Eg_count = 0.062, 0.15, 20
# y_sweep = np.linspace(Eg_start, Eg_end, Eg_count)
# y_id_str = 'Eg'
# args_to_opt= ['mu']
# arg_fix_extra = {'eta_ext':1, 'consider_nonrad':False}
# commercial_diode_ref = False
# norm_str = 'default'
# log_x = False

# Rad efficiency & Bandgap. Opt over mu.
eta_count = 80
# x_sweep = np.linspace(0.01,1, 40)
x_sweep = np.logspace(-4,0,num=eta_count,base=10)
x_id_str = 'eta_ext'
x_log = True
Eg_start = 0.02
Eg_count = 80
y_sweep = np.linspace(Eg_start, 0.3, Eg_count)
y_id_str = 'Eg'
args_to_opt = ['mu']
arg_fix_extra = {'cutoff_angle':None, 'consider_nonrad':True}
commercial_diode_ref = True
norm_str = 'log power'

alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# alg = pg.de(gen=50, ftol=1e-5)


# case_dict = {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o', 'Eg':0.094}
case_dict = {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o', 'Eg':0.094}
# {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o', 'Eg':0.1}

atm_data = atmospheric_dataset_new(case_dict['cwvstring'], location=case_dict['loc'])
case_label = case_dict['loc'] + ' ' + case_dict['cwvstring']
filename = f'PD_{case_label}_etaextlog_-4_0_{eta_count}_Egs_{Eg_start}_02_{Eg_count}.csv'
# filename = f'PD_{case_label}_cutoffangle_10_90_21_Egs_{Eg_start}_{Eg_end}_{Eg_count}.csv'

# Tc, Te = 300, 3
# atm_data = planck_law_body(Te, Ephs)
# case_label = f'T$_c$ {Tc}K, T$_e$ {Te}K'
# filename = f'PD_Te3K_Tc300K_etaextlog_-4_0_40_Egs_{Eg_start}_02_40.csv'

emitter_planck = planck_law_body(case_dict['Tskin'], atm_data.photon_energies)
combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data)

relative_to_etaext1 = False

# Generate new data
# pds_2d = []
# for y in y_sweep:
#     row = []
#     for x in x_sweep:
#     # for mu in mus:
#     #     if mu < -Eg:
#     #         pd = np.nan
#     #     else:
#     #         pd = combined_trd_env.power_density(mu, Eg, cutoff_angle=cutoff_angle)
#         arg_f = {x_id_str:x, y_id_str:y}
#         if arg_fix_extra != None:
#             arg_f.update(arg_fix_extra)
#         # best_pd = [1]
#         # while best_pd[0] > 0:
#         best_xs, best_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg=alg)
#         # print(best_pd[0])
#         row += [best_pd[0]]
#     pds_2d += [row]
# #
# pds_2d = np.asarray(pds_2d)
#
# # Save file
# np.savetxt(filename, pds_2d, delimiter = ',')

# Import existing data
pds_2d = np.loadtxt(filename, delimiter=',', dtype=float)


# Test optimizer
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


# Plot heatmap
h_ax = axs

if norm_str == 'mid0':
    # used for Eg vs V, to identify areas of TRD operation
    norm_0mid = Normalize(vmin=-10, vmax=10)
    hmap = h_ax.pcolor(x_sweep, y_sweep, pds_2d, cmap='bwr', norm=norm_0mid, shading='nearest')  # heatmap
    add_loglvls = False

elif norm_str == 'log power':
    # used for Eg vs rad efficiency, plotted logarithmically
    if relative_to_etaext1:
        # Relative to efficiency = 1 (last x)
        pds_etaext1 = pds_2d[:, -1]
        pds_relative = np.divide(pds_2d, pds_etaext1[:, np.newaxis])  # newaxis required to divide along correct axis
        C_2darray = pds_relative
        cb_label = r'P$_\mathrm{max}$ / P$_\mathrm{max}$ ($\eta_\mathrm{ext}$ = 1)'

    else:
        # Not relative - convert to + values so log scale works
        C_2darray = (-1)*np.array(pds_2d)
        cb_label = r'Power Density [W.m$^{-2}$]'

    min_a = np.min(C_2darray)
    max_a = np.max(C_2darray)
    norm_log = LogNorm(vmin=min_a, vmax=max_a)
    hmap = h_ax.pcolormesh(x_sweep, y_sweep, C_2darray, cmap='magma', shading='gouraud', norm=norm_log)
    add_loglvls = True

else:
    # default.. linear scale, no normalization.
    hmap = h_ax.pcolormesh(x_sweep, y_sweep, (-1)*pds_2d, cmap='magma', shading='gouraud')
    add_loglvls = False
    cb_label = r'Power Density [W.m$^{-2}$]'

if x_log:
    h_ax.set_xscale('log')

# Colorbar formatting
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
# cbar.ax.set_ylabel(r'Spectral Photon Flux [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]')
cbar.ax.set_ylabel(cb_label)

if add_loglvls:
    log_lvls = 10. ** np.arange(start=-7, stop=2)
    # add contour lines to heatmap
    cntr = h_ax.contour(x_sweep, y_sweep, C_2darray, levels=log_lvls, colors='black')
    # add labels to contour lines
    fmt = ticker.LogFormatterMathtext()
    fmt.create_dummy_axis()
    h_ax.clabel(cntr, cntr.levels, inline=True, fmt=fmt, fontsize=10, rightside_up=True)

    for lvl in log_lvls:
        # add ref lines in the colorbar
        cbar.ax.plot([0, 1], 2 * [lvl], '-k')


# adding commercial diodes for ref
if commercial_diode_ref:
    diodes_to_plt = [{'eta_ext':1e-2, 'Eg':0.094, 'style_args':{'marker':'o', 'color':'darkviolet', 'markersize':8}},
                     {'eta_ext':1e-1, 'Eg':0.25, 'style_args':{'marker':'o', 'mfc':'none', 'mew':2, 'color':'darkviolet', 'markersize':8}}]
    #{'eta_ext':10**(-2), 'Eg':0.175}, {'eta_ext':10**(-2), 'Eg':0.2}]
    for diode in diodes_to_plt:
        plt.plot(diode['eta_ext'], diode['Eg'], **diode['style_args'])
        ind_x = x_sweep.searchsorted(diode['eta_ext'])
        ind_y = y_sweep.searchsorted(diode['Eg'])
        cbar.ax.plot([0.5], [C_2darray[ind_y][ind_x]], **diode['style_args'])


# Labeling axes
translate_to_label = {'Eg':'Bandgap, E$_\mathrm{g}$ [eV]', 'mu':'V$_\mathrm{mpp}$ [V]',
                      'cutoff_angle': 'Cutoff Angle, $\\theta_\mathrm{c}$ [$\circ$]', 'eta_ext':'Radiative Efficiency, $\eta_\mathrm{ext}$'}
h_ax.set_xlabel(translate_to_label[x_id_str])
h_ax.set_ylabel(translate_to_label[y_id_str])
h_ax.set_xlim([x_sweep[0], x_sweep[-1]])

plt.title(f'Optimized over $\mu$, {case_label}')

plt.show()