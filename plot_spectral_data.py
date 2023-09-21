import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

cwv = 10  # column water vapor in mm, have data for 10, 24 and 54 mm

atm_dat = atmospheric_dataset(cwv)

label_to_colour = {'downwelling_0':'navy', 'downwelling_53':'lightseagreen', 'downwelling_70':'slateblue',
                   'downwelling_flux':'crimson', 'upwelling_flux':'darkmagenta', 'net_flux':'goldenrod'}

rad_units_cm = '$\mathrm{W.cm^{-2}.sr^{-1}}$'
irrad_units_cm = '$\mathrm{W.cm^{-2}}$'
rad_units_m = '$\mathrm{W.m^{-2}.sr^{-1}}$'
irrad_units_m = '$\mathrm{W.m^{-2}}$'
plot_units = [
    {'x_ret':'cm-1', 'y_ret':'W.cm-2', 'x_label':'Wavenumber, $\\tilde{v}$ [cm$^{-1}$]', 'x_array':atm_dat.wavenumbers,
     'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_cm}'+'/cm$^{-1}$]',
     'y_label_nosr':f'Spectral Irradiance, E$_e$ [{irrad_units_cm}'+'/cm$^{-1}$]'},

    {'x_ret':'um', 'y_ret':'W.m-2', 'x_label':'Wavelength, $\lambda$ [um]', 'x_array':atm_dat.wavelengths,
     'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_m}/um]',
     'y_label_nosr':f'Spectral Irradiance, E$_e$ [{irrad_units_m}/um]'},

    {'x_ret':'eV', 'y_ret':'W.m-2', 'x_label':'Photon Energy, E$_{ph}$ [eV]', 'y_label':'Spectral Radiance', 'x_array':atm_dat.photon_energies,
     'y_label_sr':f'Spectral Radiance, L$_e$ [{rad_units_m}/eV]',
     'y_label_nosr':f'Spectral Irradiance, E$_e$ [{irrad_units_m}/eV]'}
]


fig, axs = plt.subplots(2,len(plot_units), layout='tight')

radiance_labels = ['downwelling_0', 'downwelling_53', 'downwelling_70']
flux_labels = ['downwelling_flux', 'upwelling_flux', 'net_flux']

for ci, unit_set in enumerate(plot_units):
    for ri,labels in enumerate([radiance_labels, flux_labels]):
        for label in labels:
            spectral_array = atm_dat.retrieve_spectral_array(yvals=unit_set['y_ret'], xvals=unit_set['x_ret'], col_name=label)
            axs[ri][ci].plot(unit_set['x_array'],spectral_array, color=label_to_colour[label], linewidth=1.5, label=label)

    axs[0][ci].set_ylabel(unit_set['y_label_sr'])
    axs[1][ci].set_ylabel(unit_set['y_label_nosr'])
    for ri in [0,1]:
        axs[ri][ci].set_xlabel(unit_set['x_label'])


axs[0][0].legend()
axs[1][0].legend()



plt.show()



