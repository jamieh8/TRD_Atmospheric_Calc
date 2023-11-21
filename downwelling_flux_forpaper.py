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

fig, ax = plt.subplots(1)

for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str)
    Ephs = atm_data.photon_energies
    downwelling_ph_fl = atm_data.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV')
    plt.plot(Ephs, downwelling_ph_fl, c=ds['color'], label=f'{loc_str.capitalize()} {cwv_str}')

ax.set_ylabel('Spectral Photon Flux Density, F$_\mathrm{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
ax.set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
ax.legend()
# ax.set_xlim([0,0.31])
# ax.set_ylim([0,5*1e23])
ax.set_yscale('log')
add_wl_ticks(ax)


plt.show()
