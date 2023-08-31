import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

fig, axs = plt.subplots(1,1)

dat_xarray = load_file_as_xarray(cwv=10, convert_to_wavelength=True)
# in wavelength [um], [W.cm-2.sr-1/um]

dat_3D = []
for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
    dat_3D += [1e4*dat_xarray.sel(column=col_head)]  # 1e4: [cm-2] --> [m-2]

wavs = dat_xarray['wavelength']

hmap = axs.pcolor(wavs, [0,53,70], dat_3D, cmap='plasma', shading='nearest')  # heatmap

cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(r'Spectral Rad [$\mathrm{W.m^{-2}.sr^{-1}/um}$]')

axs.set_ylabel('Solid Angle [sr]')
axs.set_xlabel('Wavelength [um]')

plt.show()
