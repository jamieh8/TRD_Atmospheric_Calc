import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr

cwv = 10 # column water vapor in mm, have data fpr 10, 24 and 54 mm

filename = os.path.join('simulations_telfer', f'telfer_australia_cwv{cwv}.txt')

data = np.loadtxt(filename)

# make an xarray to add headings etc.
column_labels = [
    'wavenumber', # cm-1, centre point of bin. 0.5 cm-1 steps
    'downwelling_0', 'downwelling_53', 'downwelling_70', # RADIANCE - units of W cm-2 (cm-1)-1 sr-1
    'downwelling_flux', 'upwelling_flux', 'net_flux' # FLUX - units of W cm-2 (cm-1)-1
]
# NOTE: Helen refers to final 3 quantities as a flux but I would call this a (spectral) power density, since it has
# units of power. "Flux" for the solar spectrum is normally used to refer to the number of photons per area.

ds = xr.DataArray(data[:, 1:], coords={'wavenumber': data[:,0],
                                       'column':column_labels[1:]}, dims=['wavenumber', 'column'])

# converting wavenumber to wavelength: just 1/wavenumber. But also need to adjust flux/radiance since these
# are per cm-1 (Jacobian transformation)
ds_wl = 1e4*ds*(ds['wavenumber'])**2 # in units of (W cm-2) cm-1 -> multiply by 1e4 to get (W cm-2) um-1
ds_wl = ds_wl.rename({'wavenumber': 'wavelength'})
ds_wl.coords['wavelength'] = 1e4/ds_wl.coords['wavelength'] # convert to microns

ds_N = ds/ds['wavenumber']


fig, axs = plt.subplots(2,3, layout='tight')

label_to_colour = {'downwelling_0':'navy', 'downwelling_53':'lightseagreen', 'downwelling_70':'slateblue',
                   'downwelling_flux':'crimson', 'upwelling_flux':'darkmagenta', 'net_flux':'goldenrod'}

radiance_labels = ['downwelling_0', 'downwelling_53', 'downwelling_70']
flux_labels = ['downwelling_flux', 'upwelling_flux', 'net_flux']
for ri,labels in enumerate([radiance_labels, flux_labels]):
    for label in labels:
        axs[ri][0].plot(ds['wavenumber'], ds.sel(column=label), color=label_to_colour[label], linewidth=1.5, label=label)
        axs[ri][1].plot(ds_wl['wavelength'], ds_wl.sel(column=label), color=label_to_colour[label], linewidth=1.5, label=label)

for li, label in enumerate(['downwelling_53', 'downwelling_flux']):
    axs[li][2].plot(ds_wl['wavelength'], ds_wl.sel(column=label)*1e-4, color=label_to_colour[label], linewidth=1, label=label)


rad_units = '$\mathrm{W.cm^{-2}.sr^{-1}}$'
pd_units = '$\mathrm{W.cm^{-2}}$'
unit_xlabels = ['Wavenumber', 'Wavelength']
unit_str = ['cm$^{-1}$', 'um']
for iu in [0,1]:
    axs[0][iu].set_ylabel('Spectral Radiance [' + rad_units + f'/ {unit_str[iu]}]')
    axs[1][iu].set_ylabel(f'Spectral Power Density [{pd_units} / {unit_str[iu]}]')

    axs[0][iu].set_xlabel(unit_xlabels[iu] + f' [{unit_str[iu]}]')
    axs[1][iu].set_xlabel(unit_xlabels[iu] + f' [{unit_str[iu]}]')

axs[0][2].set_ylabel('Spectral Radiance [$\mathrm{W.m^{-2}.sr^{-1}}$ / um]')
axs[1][2].set_ylabel('Spectral Power Density [$\mathrm{W.m^{-2}}$ / um]')
axs[0][2].set_xlabel('Wavelength [um]')
axs[1][2].set_xlabel('Wavelength [um]')

axs[0][0].legend()
axs[1][0].legend()

#
# fig,ax_downwell = plt.subplots()
# ax_downwell.plot(ds_wl['wavelength'], 1e-4*ds_wl.sel(column='downwelling_flux'), label='Downwelling Flux', color=label_to_colour['downwelling_flux'])
# ax_downwell.plot(ds_wl['wavelength'], 1e-4*np.pi*ds_wl.sel(column='downwelling_53'), ':', dashes=(4, 10), label='Downwelling 53 * $\pi$', lw=1, color=label_to_colour['downwelling_53'])
# ax_downwell.set_ylabel('Spectral Power Density [ $\mathrm{W.m^{-2}}$ / um]')
# ax_downwell.set_xlabel('Wavelength [um]')
# ax_downwell.legend()


plt.show()



