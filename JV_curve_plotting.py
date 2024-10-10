import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *

set_font_opensans()

# Drawing I-V, power curves
cmap_plasma = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']
cmap_PiGrdiv = matplotlib.colormaps['PiYG']
cmap_coolwarm = matplotlib.colormaps['coolwarm']
cmap_rainbow = matplotlib.colormaps['jet']
cmap_magma = matplotlib.colormaps['magma']
hue_neg, hue_pos = 250, 15

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

Egs_planck = np.arange(0.001, 0.20, 0.02)

cutoff_angle = None

Egs_dw = np.arange(0.15, 0.20, 0.002)

to_plot = []

# Plot JV with Atmospheric Data as env
# ds = get_dataset_list()[0]
# emitter = planck_law_body(ds['Tskin'], Ephs)
# atm_dataset = atmospheric_dataset_new(ds['cwvstring'], ds['loc'], ds['Tskin'])
# combined_atmd = TRD_in_atmosphere(emitter, atm_dataset)
#
# rad_eff = 1
# Egs = np.arange(0.062, 0.2, 0.01)  # np.arange(0.155, 0.21, 0.01) #[0.094]
# Eg_opt = 0.094
# for i, Eg in enumerate(Egs):
#     # colour map around optimal Eg
#     # if Eg < Eg_opt:
#     #     coli = (Eg-Egs[0])/(Eg_opt-Egs[0]) * 0.5
#     # else:
#     #     coli = (Eg-Eg_opt)/Egs[-1] + 0.5
#     # col = cmap_coolwarm(coli)
#
#     col = cmap_plasma(i/len(Egs))
#     to_plot += [{'label':f'E$_g$ = {Eg:.3f}', 'color':col, 'linestyle':'-',
#                  'TRD in atm': combined_atmd, 'Eg':Eg,
#                  'nonrad':False, 'rad efficiency':rad_eff}]

# Plot JV with Blackbody Env
# emitter_bb = planck_law_body(300, Ephs)
# Tes = [3]
# Eg_bb = 0.1
# linestyles = ['--']
# rad_effs = [1]
# for Te, linestyle in zip(Tes, linestyles):
#     atm_bb = planck_law_body(Te, Ephs)
#     combined_atmBB = TRD_in_atmosphere(emitter_bb, atm_bb)
#     for ir, rad_eff in enumerate(rad_effs):
#         # col = cmap_purple(0.2 + 0.8 * ir / len(rad_effs))
#         col = ['black', 'red', 'blue'][ir]
#         to_plot += [{'label':f'{Te:.0f}K BB, E$_g$={Eg_bb} eV', 'color':col, 'linestyle':linestyle,
#                      'Eg':Eg_bb, 'TRD in atm': combined_atmBB,
#                      'nonrad':True, 'rad efficiency':rad_eff}]

# Plot Telfer conditions at optimal bandgaps
ds = get_dataset_list()[0:3]
for d in ds:
    dw = atmospheric_dataset_new(d['cwvstring'], d['loc'], d['Tskin'], spectral_fill_type='none', date='23dec')
    diode = planck_law_body(T=d['Tskin'])
    tia = TRD_in_atmosphere(diode,dw)
    to_plot += [{'label':'Telfer ' + d['cwvstring'], 'color':d['color'], 'linestyle':'-',
                 'TRD in atm': tia, 'Eg':d['Eg'],
                 'nonrad':False, 'rad efficiency':1}]


# Plotting options
relative_change = False  # plot change in J relative to Jsc as a function of V. True to check against Pusch et al 2019.
quadrant_focus = True  # if True, only show TRD operation quadrant
power_mult = -1

fig, axs = plt.subplots(1,2, layout='tight')
Jsc_max = 0
PD_max = 0
Voc_max = 0

# optimizer options for identifying max power point
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
    mus = np.linspace(-Eg, 0, 100)
    Js_mu = np.vectorize(TRD_in_atm_obj.current_density)(mus, Eg, cutoff_angle, eta_ext=rad_eff, consider_nonrad = nr_bool)
    pds = Js_mu*mus
    if relative_change:
        ys1 = Js_mu/Js_mu[-1]
        xs1 = q*mus/(kb * ds['Tskin'])
    else:
        xs1, ys1 = mus, Js_mu

    xs2, ys2 = mus, pds
    axs[0].plot(xs1, ys1, scheme['linestyle'], label=scheme['label'], color=i_colour)  # current
    axs[1].plot(xs2, power_mult*ys2, scheme['linestyle'], color=i_colour)  # power density

    # Test optimizer over mu
    best_xs, best_pds = get_best_pd(TRD_in_atm_obj, args_to_opt = ['mu'], args_to_fix = {'Eg':Eg, 'cutoff_angle':None, 'eta_ext':rad_eff, 'consider_nonrad':True}, alg=alg)
    axs[1].plot(best_xs['mu'], power_mult*best_pds[0], 'o', color=i_colour)

    if Js_mu[-1] > Jsc_max:
        Jsc_max = Js_mu[-1]
    if min(pds) < PD_max:
        PD_max = min(pds)

    Voc_ind = np.where(np.diff(np.signbit(Js_mu)))[0]
    if len(Voc_ind) == 0:
        Voc_max = mus[0]
    elif mus[Voc_ind[0]] < Voc_max:
        Voc_max = mus[Voc_ind[0]]


if relative_change:
    axs[0].set_xlabel('qV/(kT$_c$)')
    axs[0].set_ylabel('J/J$_{sc}$')
    axs[0].set_ylim([0, 1])
    axs[0].set_xlim([-2,0])
else:
    axs[0].set_xlabel('V [V]')
    axs[0].set_ylabel('Current Density, J [A.m$^{-2}$]')
    if quadrant_focus:
        axs[0].set_ylim([0,1.1*Jsc_max])
        # print(Voc_max)
        axs[0].set_xlim([1.1*Voc_max,0])
axs[0].legend()

axs[1].set_xlabel('V [V]')
axs[1].set_ylabel('Power Density, |JV| [W.m$^{-2}$]')

if quadrant_focus:
    axs[1].set_ylim([0,1.1*power_mult*PD_max])
    axs[1].set_xlim([1.1*Voc_max, 0])


plt.show()

