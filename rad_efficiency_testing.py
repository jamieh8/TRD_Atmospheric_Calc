import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from TRD_Atmospheric_Functions import *


# Drawing I-V, power curves
cmap_plasma = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']
cmap_purple = matplotlib.colormaps['Purples']


Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
# Egs_planck = np.arange(0.001, 0.20, 0.02)
# Egs_dw = np.arange(0.15, 0.20, 0.002)
# Egs = [0.02]
cutoff_angle = 90


to_plot = []

# cwv, Tc = 10, 296.724
Tc = 300
emitter = planck_law_body(Tc, Ephs)
rad_effs = [1,0.5,0.1]#np.arange(0.1, 1.1, 0.15)
Eg = 3*kb*Tc/q

cwv= 10
atm_dataset = atmospheric_dataset(cwv)
combined_cwv10 = TRD_in_atmosphere(emitter, atm_dataset)
for rad_eff, col in zip(rad_effs, ['dimgray', 'tomato', 'cornflowerblue']):
    to_plot += [{'label':'cwv10, $\eta_{rad}$='+f'{rad_eff:.1f}', 'color':col, 'linestyle':':',
                 'Eg':Eg, 'TRD in atm': combined_cwv10,
                 'nonrad':True, 'rad efficiency':rad_eff}]

Te_001 = Tc*0.01
Te_07 = Tc*0.7
for Te, linestyle in zip([Te_07], ['--']): #zip([Te_001, Te_07], ['-', '--']):
    atm_bb = planck_law_body(Te, Ephs)
    combined_atmBB = TRD_in_atmosphere(emitter, atm_bb)
    for ir, rad_eff in enumerate(rad_effs):
        # col = cmap_purple(0.2 + 0.8 * ir / len(rad_effs))
        col = ['black', 'red', 'blue'][ir]
        to_plot += [{'label':f'{Te:.0f} env, ' + '$\eta_{rad}$'+ f'={rad_eff:.1f}', 'color':col, 'linestyle':linestyle,
                     'Eg':Eg, 'TRD in atm': combined_atmBB,
                     'nonrad':True, 'rad efficiency':rad_eff}]


# relative_change = True
# fig, axs = plt.subplots(1,2, layout='tight')
#
# for scheme in to_plot:
#     Eg = scheme['Eg']
#     nr_bool = scheme['nonrad']
#     rad_eff = scheme['rad efficiency']
#     i_colour = scheme['color']
#     TRD_in_atm_obj = scheme['TRD in atm']
#
#     # Sweep mu for given rad efficiency
#     mus = np.linspace(-Eg, 0, 100)
#     Js_mu = np.vectorize(TRD_in_atm_obj.current_density)(mus, Eg, cutoff_angle, eta_ext=rad_eff, consider_nonrad = nr_bool)
#     pds = Js_mu*mus
#
#     if relative_change:
#         ys1 = Js_mu/Js_mu[-1]
#         xs1 = q*mus/(kb * Tc)
#     else:
#         xs1, ys1 = mus, Js_mu
#
#     xs2, ys2 = mus, pds
#     axs[0].plot(xs1, ys1, scheme['linestyle'], label = scheme['label'], color=i_colour)  # current
#     axs[1].plot(xs2, ys2, scheme['linestyle'], color=i_colour)  # power density
#
# if relative_change:
#     axs[0].set_xlabel('qV/(kT$_c$)')
#     axs[0].set_ylabel('J/J$_{sc}$')
#     axs[0].set_ylim([0, 1])
#     axs[0].set_xlim([-2,0])
# else:
#     axs[0].set_xlabel('$\mu$ [eV] / V [V]')
#     axs[0].set_ylabel('Current Density, J [A.m$^{-2}$]')
#     axs[0].set_ylim([0,2500])
# axs[0].legend()
#
# axs[1].set_xlabel('$\mu$ [eV] / V [V]')
# axs[1].set_ylabel('Power Density [W.m$^{-2}$]')
# axs[1].set_ylim([-24, 0])
#
# for ax in axs:
#     xlim_original = ax.get_xlim()
#     ax.plot(xlim_original, [0, 0], '-k', lw=0.8)
#     if relative_change == False:
#         ax.set_xlim([xlim_original[0], 0])



fig_rs, axs_rs = plt.subplots(1, 2, layout='tight')
rad_eff_sweep = np.linspace(0.01,1,10)
# rad_eff_sweep = np.logspace(-2,0,num=10,base=10)
alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)

schemes_to_plot = []
for Te, linestyle in zip([3, 210], ['-','--']):
    atm_bb = planck_law_body(Te, Ephs)
    combined_atmBB = TRD_in_atmosphere(emitter, atm_bb)
    schemes_to_plot += [{'label': f'{Te:.0f} env', 'colour': 'black', 'linestyle': linestyle,
                 'Eg': Eg, 'TRD in atm': combined_atmBB}]

for cwv, cwv_col in zip([10,24,54],['deeppink', 'steelblue', 'seagreen']):
    atm_dataset = atmospheric_dataset(cwv)
    combined_cwv = TRD_in_atmosphere(emitter, atm_dataset)
    schemes_to_plot += [{'label':f'cwv{cwv}, '+'$\eta_{rad}$=', 'colour':cwv_col, 'linestyle':'-',
                 'Eg':Eg, 'TRD in atm': combined_cwv}]

relative_change = False

for scheme in schemes_to_plot:
    TRD_in_atm_obj = scheme['TRD in atm']
    i_colour = scheme['colour']
    # Sweep rad efficiency
    Vmpps = []
    maxPds = []
    for eta in rad_eff_sweep:
        opt_xs, opt_pd = get_best_pd(TRD_in_atm_obj, args_to_opt=['mu'],
                                     args_to_fix={'Eg':Eg, 'cutoff_angle':90, 'eta_ext':eta, 'consider_nonrad':True}, alg=alg_powell)
        Vmpps += [opt_xs[0]]
        maxPds += [opt_pd[0]]

    if relative_change:
        ys1 = maxPds/maxPds[-1]
        ys2 = Vmpps/Vmpps[-1]
    else:
        ys1, ys2 = maxPds, Vmpps

    axs_rs[0].plot(rad_eff_sweep, ys1, scheme['linestyle'], label = scheme['label'], color=i_colour)  # current
    axs_rs[1].plot(rad_eff_sweep, ys2, scheme['linestyle'], label = scheme['label'], color=i_colour)  # current

for ax in axs_rs:
    ax.set_xlabel('$\eta_{ext}$')

if relative_change:
    axs_rs[0].set_ylabel('P / P($\eta$=1)')
    axs_rs[1].set_ylabel('V$_{mpp}$ / V$_{mpp}(\eta=1)$')
else:
    axs_rs[0].set_ylabel('Max Power Density [W.m$^{-2}$]')
    axs_rs[1].set_ylabel('V$_{mpp}$ [V]')


plt.show()
