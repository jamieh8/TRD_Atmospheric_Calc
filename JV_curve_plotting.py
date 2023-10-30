import matplotlib.pyplot as plt
import matplotlib
# import seaborn as sns
from TRD_Atmospheric_Functions import *

# Drawing I-V, power curves
cmap_plasma = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']
cmap_PiGrdiv = matplotlib.colormaps['PiYG']
cmap_coolwarm = matplotlib.colormaps['coolwarm']
cmap_rainbow = matplotlib.colormaps['jet']
cmap_magma = matplotlib.colormaps['magma']
hue_neg, hue_pos = 250, 15
# cmap_custom_div = sns.diverging_palette(hue_neg, hue_pos, l=70, center='dark', as_cmap=True)

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

Egs_planck = np.arange(0.001, 0.20, 0.02)

cutoff_angle = None

Egs_dw = np.arange(0.15, 0.20, 0.002)

to_plot = []

# Atmospheric Data
cwv, Tc = 10, 296.724
emitter = planck_law_body(Tc, Ephs)
atm_dataset = atmospheric_dataset(cwv)
combined_atmd = TRD_in_atmosphere(emitter, atm_dataset)

rad_eff = 1
Egs = np.arange(0.062, 0.2, 0.01) # np.arange(0.155, 0.21, 0.01)
Eg_opt = 0.094
for i, Eg in enumerate(Egs):
    # colour map around optimal Eg
    # if Eg < Eg_opt:
    #     coli = (Eg-Egs[0])/(Eg_opt-Egs[0]) * 0.5
    # else:
    #     coli = (Eg-Eg_opt)/Egs[-1] + 0.5
    # col = cmap_coolwarm(coli)
    col = cmap_plasma(i/len(Egs))
    to_plot += [{'label':f'E$_g$ = {Eg:.3f}', 'color':col, 'linestyle':'-',
                 'TRD in atm': combined_atmd, 'Eg':Eg,
                 'nonrad':False, 'rad efficiency':rad_eff}]

# Blackbody Env
# Te_001 = Tc*0.01
# Te_07 = Tc*0.7
# for Te, linestyle in zip([Te_07], ['--']): #zip([Te_001, Te_07], ['-', '--']):
#     atm_bb = planck_law_body(Te, Ephs)
#     combined_atmBB = TRD_in_atmosphere(emitter, atm_bb)
#     for ir, rad_eff in enumerate(rad_effs):
#         # col = cmap_purple(0.2 + 0.8 * ir / len(rad_effs))
#         col = ['black', 'red', 'blue'][ir]
#         to_plot += [{'label':f'{Te:.0f} env, ' + '$\eta_{rad}$'+ f'={rad_eff:.1f}', 'color':col, 'linestyle':linestyle,
#                      'Eg':Eg, 'TRD in atm': combined_atmBB,
#                      'nonrad':True, 'rad efficiency':rad_eff}]

relative_change = False
quadrant_focus = True
fig, axs = plt.subplots(1,2, layout='tight')
Jsc_max = 0
PD_max = 0
Voc_max = 0

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)
alg = alg_powell

for scheme in to_plot:
    Eg = scheme['Eg']
    nr_bool = scheme['nonrad']
    rad_eff = scheme['rad efficiency']
    i_colour = scheme['color']
    TRD_in_atm_obj = scheme['TRD in atm']

    # Sweep mu for given rad efficiency
    mus = np.linspace(-0.03*rad_eff, 0, 100)
    Js_mu = np.vectorize(TRD_in_atm_obj.current_density)(mus, Eg, cutoff_angle, eta_ext=rad_eff, consider_nonrad = nr_bool)
    pds = Js_mu*mus
    if relative_change:
        ys1 = Js_mu/Js_mu[-1]
        xs1 = q*mus/(kb * Tc)
    else:
        xs1, ys1 = mus, Js_mu

    xs2, ys2 = mus, pds
    axs[0].plot(xs1, ys1, scheme['linestyle'], label=scheme['label'], color=i_colour)  # current
    axs[1].plot(xs2, ys2, scheme['linestyle'], color=i_colour)  # power density

    # Test optimizer over mu
    best_xs, best_pds = get_best_pd(TRD_in_atm_obj, args_to_opt = ['mu'], args_to_fix = {'Eg':Eg, 'cutoff_angle':None, 'eta_ext':rad_eff, 'consider_nonrad':True}, alg=alg)
    axs[1].plot(best_xs['mu'], best_pds[0], 'o', color=i_colour)

    if Js_mu[-1] > Jsc_max:
        Jsc_max = Js_mu[-1]
    if min(pds) < PD_max:
        PD_max = min(pds)
    Voc_ind = np.where(np.diff(np.signbit(Js_mu)))[0]
    if mus[Voc_ind] < Voc_max:
        Voc_max = mus[Voc_ind]


if relative_change:
    axs[0].set_xlabel('qV/(kT$_c$)')
    axs[0].set_ylabel('J/J$_{sc}$')
    axs[0].set_ylim([0, 1])
    axs[0].set_xlim([-2,0])
else:
    axs[0].set_xlabel('$\mu$ [eV] / V [V]')
    axs[0].set_ylabel('Current Density, J [A.m$^{-2}$]')
    if quadrant_focus:
        axs[0].set_ylim([0,1.1*Jsc_max])
        axs[0].set_xlim([1.1*Voc_max,0])
axs[0].legend()

axs[1].set_xlabel('$\mu$ [eV] / V [V]')
axs[1].set_ylabel('Power Density [W.m$^{-2}$]')
if quadrant_focus:
    axs[1].set_ylim([1.1*PD_max,0])
    axs[1].set_xlim([1.1 * Voc_max, 0])




plt.show()