from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

datasets = [
    {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o'},
    {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o'},
    {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o'},

    # {'loc':'california', 'cwvstring':'low', 'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'mid', 'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'high', 'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'symbol': 's'},

    # {'loc':'tamanrasset', 'cwvstring':'low', 'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'mid', 'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'high', 'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'symbol':'^'}
    ]

fig = plt.figure()
gs1 = fig.add_gridspec(2, 3, width_ratios=[2,2,3], wspace=0.5, hspace=0.5)
# gs2 = fig.add_gridspec(1, 3, width_ratios=[5,1,2], wspace=0.05)

ax_latmap = fig.add_subplot(gs1[:,0])
ax_hist = fig.add_subplot(gs1[0,1])
ax_scatter = fig.add_subplot(gs1[1,1])
ax_dwf = fig.add_subplot(gs1[:,2])
axs = [ax_latmap, ax_hist, ax_scatter, ax_dwf]

# Histogram and scatter plot
copied_data = np.loadtxt('histogram-plot-data.csv', skiprows=1, delimiter=',')
edges = np.arange(0,copied_data[-1,0]+1,1)
ax_hist.stairs(copied_data[:,1], edges=edges, fill=True, fc='silver')

vert_line_dicts = [
    {'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'ls': 'dotted', 'symbol': 's'},

    {'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'ls':'solid', 'symbol':'o'},
    {'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','ls':'solid', 'symbol':'o'},
    {'tcwv':70.51, 'Tskin':299.86, 'color':'teal','ls':'solid', 'symbol':'o'},

    {'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'ls':'dashed', 'symbol':'^'},
    {'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'ls':'dashed', 'symbol':'^'},
    {'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'ls':'dashed', 'symbol':'^'}]

for line_dict in vert_line_dicts:
    tcwv = line_dict['tcwv']
    Tskin = line_dict['Tskin']
    insert_index = np.searchsorted(edges, line_dict['tcwv'])
    try:
        yval = copied_data[:,1][insert_index-1]
    except:
        yval = 0
    # axs[0].plot(2*[line_dict['tcwv']], [0,4], ls=line_dict['ls'], c=line_dict['color'])
    ax_hist.annotate('', xy=(tcwv, yval), xytext=(tcwv,yval+0.5), arrowprops = dict(arrowstyle='->', color='black'))
    ax_hist.plot([tcwv], [yval+0.55], line_dict['symbol'], c=line_dict['color'])

    ax_scatter.plot([tcwv], [Tskin], line_dict['symbol'], c=line_dict['color'], markersize=7)

ax_scatter.set_xlabel('TCWV [mm]')
ax_scatter.set_ylabel('Skin Temp [K]')
ax_scatter.minorticks_on()
ax_scatter.grid()

ax_hist.set_xlabel('TCWV [mm]')
ax_hist.set_ylabel('Occurence [%]')
ax_hist.set_ylim([0,4])



# Downwelling flux
for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str)
    Ephs = atm_data.photon_energies
    downwelling_ph_fl = atm_data.retrieve_spectral_array(yvals='W.cm-2', xvals='cm-1')
    ax_dwf.plot(Ephs, downwelling_ph_fl, c=ds['color'], label=f'{loc_str.capitalize()} {cwv_str}')

ax_dwf.set_ylabel('Spectral Photon Flux Density, F$_\mathrm{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
# ax.set_ylabel('Spectral Irradiance, F$_e$ [W.cm$^{-2}$/cm$^{-1}$]')
ax_dwf.set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
ax_dwf.legend(loc='lower left')
# ax.set_xlim([0,0.31])
# ax.set_ylim([0,5*1e23])
ax_dwf.set_yscale('log')

add_wl_ticks(ax_dwf)
# add_wn_ticks(ax)

for letter, ax in zip(['a','b','c','d'], axs):
    ax.set_title(letter+')', loc='left', fontsize=15, fontweight='bold', y=1.05, x=-0.1)

plt.show()
