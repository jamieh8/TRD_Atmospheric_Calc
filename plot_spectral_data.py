import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *


ds_lsts = [get_dataset_list()[0:3], get_dataset_list()[3:6], get_dataset_list()[6:9]]


# label_to_colour = {'downwelling_flux':'darkorange', 'upwelling_flux':'darkmagenta', 'net_flux':'goldenrod'}

fig, axs = plt.subplots(2,3, layout='tight')
for ci, ds_lst in enumerate(ds_lsts):
    for ds in ds_lst:
        atm_dat = atmospheric_dataset_new(cwv=ds['cwvstring'], location=ds['loc'], Tskin=ds['Tskin'], date='23dec')
        ys = atm_dat.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV', col_name='downwelling_flux')
        axs[0][ci].plot(atm_dat.photon_energies, ys, c=ds['color'], label=ds['loc']+' '+ds['cwvstring'], alpha=0.8)

        ys = atm_dat.retrieve_spectral_array(yvals='W.m-2', xvals='um', col_name='downwelling_flux')
        axs[1][ci].plot(atm_dat.wavelengths, ys, c=ds['color'], alpha=0.8)

for ax in axs[0][:]:
    ax.set_yscale('log')
    ax.set_xlabel('Photon energy, $E_g$ [eV]')
    ax.set_ylabel('Spectral Photon Flux Density, F$_{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
    ax.legend()
    ax.grid()

for ax in axs[1][:]:
    ax.set_xlim([0,30])
    ax.set_xlabel('Wavelength [um]')
    ax.set_ylabel('Spectral Irradiance, F$_e$ [W.m$^{-2}$/um]')
    ax.grid()


# rad_units_cm = '$\mathrm{W.cm^{-2}.sr^{-1}}$'
# irrad_units_cm = '$\mathrm{W.cm^{-2}}$'
# rad_units_m = '$\mathrm{W.m^{-2}.sr^{-1}}$'
# irrad_units_m = '$\mathrm{W.m^{-2}}$'
# plot_units = [
#     {'x_ret':'cm-1', 'y_ret':'W.cm-2', 'x_label':'Wavenumber, $\\tilde{v}$ [cm$^{-1}$]', 'x_array':atm_dat.wavenumbers,
#      'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_cm}'+'/cm$^{-1}$]',
#      'y_label_nosr':f'Spectral Irradiance, F$_e$ [{irrad_units_cm}'+'/cm$^{-1}$]'},
#
#     {'x_ret':'um', 'y_ret':'W.m-2', 'x_label':'Wavelength, $\lambda$ [um]', 'x_array':atm_dat.wavelengths,
#      'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_m}/um]',
#      'y_label_nosr':f'Spectral Irradiance, F$_e$ [{irrad_units_m}/um]'},
#
#     {'x_ret':'eV', 'y_ret':'W.m-2', 'x_label':'Photon Energy, E$_{ph}$ [eV]', 'y_label':'Spectral Radiance', 'x_array':atm_dat.photon_energies,
#      'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_m}/eV]',
#      'y_label_nosr':f'Spectral Irradiance, F$_e$ [{irrad_units_m}/eV]'}
# ]

# fig, axs = plt.subplots(2,len(plot_units), layout='tight')
# flux_labels = ['downwelling_flux']#, 'upwelling_flux', 'net_flux']
# radiance_labels = [f'downwelling_{theta}' for theta in zenith_angles]

# cmap = matplotlib.colormaps['turbo']
# for i, rad_lab in enumerate(radiance_labels):
#     label_to_colour.update({rad_lab:cmap(i/(len(zenith_angles)-1))})
# for ci, unit_set in enumerate(plot_units):
#     for ri,labels in enumerate([radiance_labels, flux_labels]):
#         for label in labels:
#             spectral_array = atm_dat.retrieve_spectral_array(yvals=unit_set['y_ret'], xvals=unit_set['x_ret'], col_name=label)
#             axs[ri][ci].plot(unit_set['x_array'],spectral_array, color=label_to_colour[label], linewidth=1.5, label=label)
#
#     axs[0][ci].set_ylabel(unit_set['y_label_sr'])
#     axs[1][ci].set_ylabel(unit_set['y_label_nosr'])
#     for ri in [0,1]:
#         axs[ri][ci].set_xlabel(unit_set['x_label'])
# axs[0][0].legend()
# axs[1][0].legend()



plt.show()



