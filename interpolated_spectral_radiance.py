import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

fig, axs = plt.subplots(1,2, layout='tight')

dat_xarray = load_file_as_xarray(cwv=10, convert_to_wavelength=True)
# in wavelength [um], [W.cm-2.sr-1/um]

dat_3D = []
for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
    dat_3D += [1e4*dat_xarray.sel(column=col_head)]  # 1e4: [cm-2] --> [m-2]
dat_3D = np.array(dat_3D)
wavs = dat_xarray['wavelength']


# interpolate data at different angles
angle_array = np.linspace(0,90,10)
dat_3D_int = []
for wi in range(len(wavs)):
    # for each wavelength, interpolate values for angles specified
    dat_3D_int += [np.interp(x=angle_array, xp=[0,53,70], fp=dat_3D[:,wi])]

dat_3D_int = np.array(dat_3D_int)  # rows correspond to angles
dat_3D_int = dat_3D_int.transpose()  # rows correspond to wavelengths

# make heatmap
h_ax = axs[0]
hmap = h_ax.pcolor(wavs, angle_array, dat_3D_int, cmap='plasma', shading='nearest')  # heatmap
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(r'Spectral Rad [$\mathrm{W.m^{-2}.sr^{-1}/um}$]')
h_ax.set_ylabel('Solid Angle [sr]')
h_ax.set_xlabel('Wavelength [um]')


# take slices of heatmap, plot separately as line plot
slice_ax = axs[1]
slice_wavelengths = [6,9.5,13,18]
s_idxs = len(wavs) - np.searchsorted(np.flip(wavs),slice_wavelengths)  # indexes of slice wavelengths
for si, swav in zip(s_idxs, slice_wavelengths):
    spectral_rad_slice = dat_3D_int[:,si]
    slice_ax.plot(angle_array, spectral_rad_slice, label=f'{swav} um')

    h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--')

slice_ax.set_xlabel('Solid Angle [sr]')
slice_ax.set_ylabel('Spectral Rad [$\mathrm{W.m^{-2}.sr^{-1}/um}$]')
slice_ax.legend()

plt.show()
