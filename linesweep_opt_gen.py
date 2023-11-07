import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

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
# args_to_opt = ['mu', 'cutoff_angle']
# x_vals_sweep = np.linspace(0.062,0.2,20)  # bandgaps, in [eV]
# arg_sweep_str = 'Eg'
# arg_fix_extra = {'eta_ext':1, 'consider_nonrad':False}
# xlabel_str = 'Bandgap, E$_g$ [eV]'

# Sweep cutoff angle, optimize over mu (with Eg defined per case)
args_to_opt = ['mu']
x_vals_sweep = np.arange(10,91,5)
arg_sweep_str = 'cutoff_angle'
xlabel_str = 'Cutoff Angle, $\\theta_c$ [$\circ$]'
fix_Eg = True
arg_fix_extra = {}

to_plot = []


# Setup cases to plot
# old_file_dicts = [
#     {'cwv':10, 'Tc':296.724, 'col':'deeppink', 'Eg':0.094, 'new':False},
#     {'cwv':24, 'Tc':304.868, 'col':'steelblue', 'Eg':0.094, 'new':False},
#     {'cwv':54, 'Tc':303.512, 'col':'seagreen', 'Eg':0.1, 'new':False}
# ]

new_file_dicts = [
    # {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o', 'Eg':0.094},
    {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o', 'Eg':0.094},
    # {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o', 'Eg':0.1},

    # {'loc':'california', 'cwvstring':'low', 'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'mid', 'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'high', 'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'symbol': 's'},
    #
    # {'loc':'tamanrasset', 'cwvstring':'low', 'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'mid', 'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'high', 'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'symbol':'^'}
    ]

for AD_dict in new_file_dicts:
    try:
        Eg = AD_dict['Eg']
    except:
        Eg = 0.094
    cwv = AD_dict['tcwv']
    cwv_str = AD_dict['cwvstring']
    loc_str = AD_dict['loc']
    atm_data = atmospheric_dataset_new(cwv_str, loc_str)
    TRD = planck_law_body(T=AD_dict['Tskin'], Ephs=atm_data.photon_energies)
    comb_TRDenv = TRD_in_atmosphere(TRD, atm_data)
    if fix_Eg:
        arg_fix_extra.update({'Eg':Eg})  #<--- comment out if Eg is to be optimized
    to_plot += [{'TRD in env':comb_TRDenv, 'label':f'{loc_str} {cwv_str}', 'colour':AD_dict['color'], 'arg_fix_extra': arg_fix_extra}]




# TRD_300K = planck_law_body(T=300, Ephs=Ephs)
# for Te, colT in zip([3,200],['black','dimgrey']):
#     env_T = planck_law_body(T=Te, Ephs=Ephs)
#     comb_TRDenv = TRD_in_atmosphere(TRD_300K, env_T)
#     to_plot += [{'TRD in env':comb_TRDenv, 'label':f'300K/{Te}K, E$_g$=0.02', 'colour':colT, 'arg_fix_extra': {'Eg':0.02}}]


# Define subplots based on optimization
fig, axs = plt.subplots(1,len(args_to_opt)+1, layout='tight')
relative_change = True

for case in to_plot:
    save_array = []
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
        opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=args_to_opt, args_to_fix=arg_f, alg = alg_powell)
        max_pds += [opt_pd[0]]
        opt_vals += [[opt_xs.get(str_lbl) for str_lbl in args_to_opt]]  # ensures order of values in opt_vals matches order of args_to_opt.!

        # print(f'optimized xs: {opt_xs}')
        # print(f'fixed: {arg_f}')
        print(f'max PD: {opt_pd[0]}')
        save_array += [[xval, opt_pd[0]]]
        # optres_90, pd_90 = get_best_pd(combined_trd_env, args_to_opt=['mu'], args_to_fix={arg_sweep_str: xval, 'cutoff_angle':None}, alg = alg_powell)
        # max_pds_90 += [pd_90[0]]
        # mus_90 += [optres_90['mu']]

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
    # axs[args_to_opt.index('cutoff_angle')+1].plot(x_vals_sweep, len(x_vals_sweep)*[90],'--', c=case['colour'])
    # axs[args_to_opt.index('mu')+1].plot(x_vals_sweep, mus_90, '--', c=case['colour'])
    lbl = case['label']
    filename = f'{lbl}_sweep_{arg_sweep_str}.csv'
    # np.savetxt(fname = filename, X = np.array(save_array), delimiter=',')


pd_ax = axs[0]
if relative_change:
    # pd_ax.set_ylabel(f'Relative Power Density, P/P({x_vals_sweep[-1]:.0f})')
    pd_ax.set_ylabel('Relative Power Density, P($\\theta_c$)/P($90^\circ$)')
else:
    pd_ax.set_ylabel(r'Power Density [W.m$^{-2}$]')
pd_ax.set_xlabel(xlabel_str)
pd_ax.legend()

translate_to_label = {'Eg':'Bandgap, E$_g$ [eV]', 'mu':'V$_{mpp}$ [V]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}
for opi, opt_label in enumerate(args_to_opt):
    axs[opi+1].set_ylabel(translate_to_label[opt_label])
    axs[opi+1].set_xlabel(xlabel_str)

axs[0].plot([0,90],[1,1], '--k', lw=1)
axs[0].set_xlim([10,90])

plt.show()


