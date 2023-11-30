import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, text
import matplotlib.gridspec as gridspec

# Heatmap
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


# Import existing data
pds_2d = np.loadtxt(filename, delimiter=',', dtype=float)



# Radiative efficiency testing:
comparison_lst = []
args_to_opt = ['mu']

datasets = get_dataset_list()

for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str)
    Ephs = atm_data.photon_energies
    line_format_dct = {'color': ds['color'], 'linestyle': 'solid'}
    scatter_format = {'c': ds['color'], 'marker': ds['symbol'], 'markersize':8}
    emitter_planck = planck_law_body(T=ds['Tskin'], Ephs=Ephs)
    comparison_lst += [{'label': f'{loc_str} {cwv_str}', 'colour':ds['color'], 'scatter format':scatter_format,
                        'cwv':ds['tcwv'],
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, atm_data)}]



emitter_planck = planck_law_body(T=300)
env_planck = planck_law_body(T=3)
comparison_lst += [{'label': f'3K', 'colour':'black','scatter format': {'c':'k', 'marker':'o', 'markersize':8},
                        'cwv':0,
                        'TRD in atm':TRD_in_atmosphere(emitter_planck, env_planck)}]

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)



# fig, axs = plt.subplots(1,3,layout='tight', sharey=False, width_ratios=[4,1,2])
fig = plt.figure()

gs1 = fig.add_gridspec(1, 3, width_ratios=[5,1,2], wspace=0.5)
gs2 = fig.add_gridspec(1, 3, width_ratios=[5,1,2], wspace=0.05)

ax1 = fig.add_subplot(gs1[0])
ax2 = fig.add_subplot(gs2[1])
ax3 = fig.add_subplot(gs2[2], sharey = ax2)
axs = [ax1,ax2,ax3]

# Plot heatmap
h_ax = axs[0]

# used for Eg vs rad efficiency, plotted logarithmically
C_2darray = (-1)*np.array(pds_2d)
cb_label = r'Power Density [W.m$^{-2}$]'

min_a = 1e-7#np.min(C_2darray)
max_a = 3#np.max(C_2darray)
norm_log = LogNorm(vmin=min_a, vmax=max_a)
hmap = h_ax.pcolormesh(x_sweep, y_sweep, C_2darray, cmap='magma', shading='gouraud', norm=norm_log)
add_loglvls = True

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
    diodes_to_plt = [{'eta_ext':1, 'Eg':0.094, 'style_args':{'marker':'o', 'color':'darkviolet', 'markersize':8}},
                      {'eta_ext':1e-2, 'Eg':0.094, 'style_args':{'marker':'o', 'color':'darkviolet', 'markersize':8, 'fillstyle':'right'}},
                     {'eta_ext':1e-1, 'Eg':0.25, 'style_args':{'marker':'o', 'mfc':'none', 'mew':1.5, 'color':'darkviolet', 'markersize':8}}]
    #{'eta_ext':10**(-2), 'Eg':0.175}, {'eta_ext':10**(-2), 'Eg':0.2}]
    for diode in diodes_to_plt:
        h_ax.plot(diode['eta_ext'], diode['Eg'], 'o', markersize=12, c='white')
        h_ax.plot(diode['eta_ext'], diode['Eg'], **diode['style_args'])
        ind_x = x_sweep.searchsorted(diode['eta_ext'])
        ind_y = y_sweep.searchsorted(diode['Eg'])
        cbar.ax.plot([0.5], [C_2darray[ind_y][ind_x]], 'o', markersize=12, c='white')
        cbar.ax.plot([0.5], [C_2darray[ind_y][ind_x]], **diode['style_args'])



h_ax.set_xlabel('Radiative Efficiency, $\eta_\mathrm{ext}$')
h_ax.set_ylabel('Bandgap, E$_\mathrm{g}$ [eV]')
h_ax.set_xlim([x_sweep[0], x_sweep[-1]])


diodes = [{'Eg':0.094, 'eta':1, 'styleargs':{}}, {'Eg':0.094, 'eta':0.01, 'styleargs':{'fillstyle':'right'}}, {'Eg':0.25, 'eta':0.1, 'styleargs':{'mfc':'none', 'mew':1.5}}]
ax_scatter = axs[1]
for diode in diodes:
    for case in comparison_lst:
        opt_xs, opt_pd = get_best_pd(case['TRD in atm'],
                                     args_to_opt=args_to_opt,
                                     args_to_fix={'cutoff_angle':None, 'eta_ext':diode['eta'], 'Eg':diode['Eg'], 'consider_nonrad':True},
                                     alg=alg_powell)

        style_args = diode['styleargs']
        style_args.update(case['scatter format'])
        ax_scatter.plot(case['cwv'], opt_pd[0]*(-1), **style_args)


ax_scatter.set_yscale('log')
ax_scatter.set_xlabel('TCWV [mm]')
ax_scatter.set_ylabel('Power Density [W.m$^{-2}$]')
ax_scatter.minorticks_on()

reflines = [{'y':29, 'string':'yearly avg solar', 'xl':-5, 'xr':10, 'xt':10, 'arrow args':{'arrowstyle':'<-', 'color':'orangered'}, 'text args':{'color':'orangered', 'va':'center', 'ha':'left'}},
            {'y':51, 'string':'300K to 3K limit', 'xl':-5, 'xr':10, 'xt':10, 'arrow args':{'arrowstyle':'<-', 'color':'black'}, 'text args':{'color':'black', 'va':'center', 'ha':'left'}}]
for rl in reflines:
    ax_scatter.annotate(text='', xy=(rl['xr'],rl['y']), xytext=(rl['xl'],rl['y']), arrowprops=rl['arrow args'])
    ax_scatter.plot([0,80],2*[rl['y']], color='white', lw=1)
    ax_scatter.text(s=rl['string'], x=rl['xt'], y=rl['y'], **rl['text args'])



# Plot power magnitude examples
sample_area = 10  # [m-2]
power_magnitudes_guides = [{'label':'Average single person household \nin Australia for 1 day', 'PD':8*1e3 / (sample_area*24)},
                           {'label':'800W microwave for 15 mins', 'PD':800*0.25 / (sample_area*24)},
                           {'label':'10W LED bulb for 5h', 'PD':10*5 / (sample_area*24)},
                           {'label':'5W phone charger for 2h', 'PD':10 / (sample_area*24)},
                           {'label':'60W fridge for 24h', 'PD':60*24 / (sample_area*24)},]
ax_powerguide = axs[2]
for pguide in power_magnitudes_guides:
    ax_powerguide.annotate(pguide['label'], xy=(0,pguide['PD']), xytext=(0.1, pguide['PD']), arrowprops=dict(arrowstyle="->"), va='center')
    # ax_powerguide.text(s=pguide['label'], x=0.05, y=pguide['PD'], ha='left', va='bottom')
    # print(pguide['label'] + ' ' + str(pguide['PD']*24))
ax_powerguide.set_title(f'   Over 24h, {sample_area} m$^2$ can power:', loc='left')
ax_powerguide.axis('off')

# Add legend in second plot
custom_muopt_legend = []
for diode in diodes:
    Eg, eta = diode['Eg'], diode['eta']
    custom_muopt_legend += [Line2D([0],[0], linestyle='none', **diode['styleargs'],
                                   label='E$_\mathrm{g}$=' + f'{Eg} eV' + ', $\eta_\mathrm{ext}=$' + f'{eta*100:.0f}%')]

ax_powerguide.legend(handles=custom_muopt_legend, loc='lower left')

ax_scatter.grid(axis='y')
ax_scatter.set_xlim([-5,75])
ax_scatter.set_ylim([1e-4, 1e2])

# text.Annotation('figure points text', xy=(0.5,0.5), xycoords='figure fraction',fontsize=15, fontweight='bold')
ax_scatter.set_title('b)', loc='left', fontsize=15, fontweight='bold', pad=15, x=-0.3)
h_ax.set_title('a)', loc='left', fontsize=15, fontweight='bold', pad=15, x=-0.15)

plt.show()