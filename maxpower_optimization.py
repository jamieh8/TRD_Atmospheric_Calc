import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar
from matplotlib.colors import Normalize


# ----------------- Line/Sweep Optimization Plotter -----------------

# Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
#
# # cwv, T = 10, 296.724
# # cwv, T = 24, 304.868
# # cwv, T = 54, 303.512
#
# to_plot = []
# to_plot += [{'cwv':10, 'Tc':296.724, 'label':'cwv10', 'colour':'deeppink'}]
# to_plot += [{'cwv':24, 'Tc':304.868, 'label':'cwv24', 'colour':'steelblue'}]
# to_plot += [{'cwv':54, 'Tc':303.512, 'label':'cwv54', 'colour':'seagreen'}]
#
# alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# # alg = pg.de(gen=50, ftol=1e-5)
#
# # Sweep cutoff angle, optimize over Eg and mu
# # args_to_opt = ['Eg', 'mu']
# # x_vals_sweep = np.arange(10,95,10)
# # arg_sweep_str = 'cutoff_angle'
# # arg_fix_extra = {}
# # xlabel_str = 'Cutoff Angle [$\circ$]'
#
# # Sweep Eg, fix cutoff angle 90, optimize over mu
# # args_to_opt = ['mu']
# # x_vals_sweep = np.linspace(0.062,0.15,20)  # bandgaps, in [eV]
# # arg_sweep_str = 'Eg'
# # arg_fix_extra = {'cutoff_angle':90}
# # xlabel_str = 'Bandgap, E$_g$ [eV]'
#
# # Sweep Eg, optimize over mu and cutoff angle
# args_to_opt = ['mu', 'cutoff_angle']
# x_vals_sweep = np.linspace(0.062,0.15,20)  # bandgaps, in [eV]
# arg_sweep_str = 'Eg'
# arg_fix_extra = {}
# xlabel_str = 'Bandgap, E$_g$ [eV]'
#
#
# # Define subplots based on optimization
# fig, axs = plt.subplots(1,len(args_to_opt)+1, layout='tight')
#
# for case in to_plot:
#     atm_data = atmospheric_dataset(case['cwv'])
#     emitter_planck = planck_law_body(T=case['Tc'])
#     combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data, Ephs)
#
#     max_pds = []
#     opt_vals = []
#     for xval in x_vals_sweep:
#         arg_f = {arg_sweep_str: xval}
#         if arg_fix_extra != None:
#             arg_f.update(arg_fix_extra)
#         opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg = alg)
#         max_pds += [opt_pd]
#         opt_vals += [opt_xs]
#
#     axs[0].plot(x_vals_sweep, max_pds, c=case['colour'], label=case['label'])
#
#     # Secondary plot showing x values corresponding to best P.D.
#     opt_vals_t = np.transpose(np.array(opt_vals))
#     for ir, row in enumerate(opt_vals_t):
#         axs[ir+1].plot(x_vals_sweep, row, c=case['colour'], label=case['label'])
#
# pd_ax = axs[0]
# pd_ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
# pd_ax.set_xlabel(xlabel_str)
# pd_ax.legend()
#
# translate_to_label = {'Eg':'Bandgap, E$_g$ [eV]', 'mu':'V [V] / $\mu$ [eV]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}
# for opi, opt_label in enumerate(args_to_opt):
#     axs[opi+1].set_ylabel(translate_to_label[opt_label])
#     axs[opi+1].set_xlabel(xlabel_str)



# ----------------- Heatmap for testing optimizer -----------------
fig, axs = plt.subplots(1,1)
Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
Egs = np.linspace(0.062, 0.2, 100)
mus = np.linspace(-0.1, 0, 100)
# cutoff_angle = 70
cutoff_angles = np.linspace(10,90,20)

# filename = f'PD_cutoff{cutoff_angle}_Egs_0.062_02_100_mus_-01_0_100.csv'
# filename = f'PD_optmu_Egs_0.062_02_100_cutoffangle_'

cwv, T = 10, 296.724
atm_data = atmospheric_dataset(cwv)
emitter_planck = planck_law_body(T)
combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data, Ephs)

# # Generate new data
pds_2d = []
for Eg in Egs:
    row = []
    for mu in mus:
        if mu < -Eg:
            pd = np.nan
        else:
            pd = combined_trd_env.power_density(mu, Eg, cutoff_angle=cutoff_angle)
        row += [pd]
    pds_2d += [row]
# np.savetxt(filename, np.asarray(pds_2d), delimiter = ',')

# Import existing data
# pds_2d = np.loadtxt(filename, delimiter=',', dtype=float)

# Test optimizer
# alg = pg.scipy_optimize(method='Powell', tol=1e-5)
# # alg = pg.de(gen=50, ftol=1e-5)
# Egs_opt = []
# Vs_opt = []
# pds_opt = []
# for i in range(5):
#     opt_res = get_best_pd(combined_trd_env, args_to_opt={'Eg':0, 'mu':0}, args_to_fix={'cutoff_angle':cutoff_angle}, alg = alg)
#     plt.plot([opt_res[0][1]], [opt_res[0][0]], 'x', c='black')
#
#     Egs_opt += [opt_res[0][1]]
#     Vs_opt += [opt_res[0][0]]
#     pds_opt += [opt_res[1][0]]
#
# for str, array in zip(['Egs','Vs','PDs'], [Egs_opt, Vs_opt, pds_opt]):
#     print(str)
#     for elem in array:
#         print(elem)

# Plot heatmap
norm_0mid = Normalize(vmin=-10, vmax=10)
h_ax = axs
hmap = h_ax.pcolor(mus, Egs, pds_2d, cmap='bwr', norm=norm_0mid, shading='nearest')  # heatmap
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
# cbar.ax.set_ylabel(r'Spectral Photon Flux [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]')
cbar.ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
h_ax.set_xlabel('V [V] / $\mu$ [eV]')
h_ax.set_ylabel('Bandgap, E$_g$ [eV]')


plt.show()