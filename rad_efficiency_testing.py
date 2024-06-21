import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.lines import Line2D

# Plotting power output vs radiative efficiency (fixed bandgap). Used to check against Pusch et al 2019.

# ------------------------ Sweep eta_ext ------------------------
Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
cutoff_angle = None


schemes_to_plot = []

# Blackbody environment
for Te, linestyle in zip([3, 210], ['-','--']):
    atm_bb = planck_law_body(Te, Ephs)
    emitter = planck_law_body(T=300)
    combined_atmBB = TRD_in_atmosphere(emitter, atm_bb)
    schemes_to_plot += [{'label': f'{Te:.0f} env', 'colour': 'black', 'linestyle': linestyle,
                 'Eg': 0.1, 'TRD in atm': combined_atmBB}]


# LBLRTM modelling environment
dss = [get_dataset_list()[0]]

for ds in dss:
    label = ds['loc'] + ' ' + ds['cwvstring']
    atm_dataset = atmospheric_dataset_new(cwv=ds['cwvstring'], location=ds['loc'], Tskin=ds['Tskin'], date='23dec')
    emitter = planck_law_body(T=ds['Tskin'])
    combined_cwv = TRD_in_atmosphere(emitter, atm_dataset)
    schemes_to_plot += [{'label':label,'colour':ds['color'], 'linestyle':'-', 'Eg':ds['Eg'], 'TRD in atm': combined_cwv}]


log_eta = True
if log_eta:
    rad_eff_sweep = np.logspace(-4, 0, num=10, base=10)
else:
    rad_eff_sweep = np.linspace(0.01,1,10)

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
args_to_opt = ['mu']


fig_rs, axs_rs = plt.subplots(1, 1+len(args_to_opt), layout='tight')
relative_change = True

for scheme in schemes_to_plot:
    TRD_in_atm_obj = scheme['TRD in atm']
    i_colour = scheme['colour']

    # Sweep rad efficiency
    maxPds = []
    opt_vals = []
    for eta in rad_eff_sweep:
        opt_xs, opt_pd = get_best_pd(TRD_in_atm_obj, args_to_opt=args_to_opt,
                                     args_to_fix={'cutoff_angle':90, 'eta_ext':eta, 'consider_nonrad':True, 'Eg':scheme['Eg']},
                                     alg=alg_powell)
        print(opt_xs)
        opt_vals += [list(opt_xs.values())]
        maxPds += [opt_pd[0]]

    if relative_change:
        ys1 = maxPds/maxPds[-1]
        # ys2 = Vmpps/Vmpps[-1]
    else:
        ys1 = maxPds

    if log_eta:
        axs_rs[0].loglog(rad_eff_sweep, ys1, scheme['linestyle'], label = scheme['label'], color=i_colour)  # current
    else:
        axs_rs[0].plot(rad_eff_sweep, ys1, scheme['linestyle'], label=scheme['label'], color=i_colour)  # current

    # Secondary plot showing x values corresponding to best P.D.
    opt_vals_t = np.transpose(np.array(opt_vals))
    for ir, row in enumerate(opt_vals_t):
        if log_eta:
            axs_rs[ir+1].loglog(rad_eff_sweep, np.abs(row), scheme['linestyle'], c=scheme['colour'], label=scheme['label'])
        else:
            axs_rs[ir + 1].loglog(rad_eff_sweep, row, scheme['linestyle'], c=scheme['colour'],label=scheme['label'])


for ax in axs_rs:
    ax.set_xlabel('$\eta_{ext}$')

translate_to_label = {'Eg':'Bandgap, E$_\mathrm{g}$ [eV]', 'mu':'V$_\mathrm{mpp}$ [V]', 'cutoff_angle': 'Cutoff Angle [$\circ$]'}

if 'mu' in args_to_opt:
    args_to_opt.remove('mu')
    args_to_opt += ['mu']
for opi, opt_label in enumerate(args_to_opt):
    axs_rs[opi+1].set_ylabel(translate_to_label[opt_label])

if relative_change:
    axs_rs[0].set_ylabel('P / P($\eta$=1)')
else:
    axs_rs[0].set_ylabel('Max Power Density [W.m$^{-2}$]')


plt.show()
