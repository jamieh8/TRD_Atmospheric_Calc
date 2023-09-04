import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

# make heatmap
# xstr = 'photon energies'
# int_dict = interpolate_spectralrad_by_angle(cwv=10, x_vals=xstr, number_of_angles=100)
# xs = int_dict[f'{xstr} (columns)']
# dat_3D_int = int_dict['int rad']
# angle_array = int_dict['angles (rows)']
# fig, axs = plt.subplots(1,2, layout='tight')
#
# h_ax = axs[0]
# hmap = h_ax.pcolor(xs, angle_array, dat_3D_int, cmap='plasma', shading='nearest')  # heatmap
# cbar = plt.colorbar(hmap)
# cbar.ax.tick_params(labelsize=10)
# cbar.ax.set_ylabel(r'Spectral Rad [$\mathrm{W.m^{-2}.sr^{-1}/um}$]')
# h_ax.set_ylabel('Solid Angle [sr]')
# h_ax.set_xlabel(xstr)
#
#
# # take slices of heatmap, plot separately as line plot
# slice_ax = axs[1]
# if xstr == 'wavelengths':
#     slice_xs = [6,9.5,13,18]
#     s_idxs = len(xs) - np.searchsorted(np.flip(xs), slice_xs)  # indexes of slice wavelengths
# elif xstr == 'photon energies':
#     slice_xs = [0.08, 0.12, 0.15, 0.2]
#     s_idxs = np.searchsorted(xs, slice_xs)
#
# for si, swav in zip(s_idxs, slice_xs):
#     # plot spectral rad vs angle for a particular photon energy / wavelength
#     spectral_rad_slice = dat_3D_int[:,si]
#     slice_ax.plot(angle_array, spectral_rad_slice, label=f'{swav}')
#
#     h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--')
#
# slice_ax.set_xlabel('Solid Angle [sr]')
# slice_ax.set_ylabel('Spectral Rad [$\mathrm{W.m^{-2}.sr^{-1}/um}$]')
# slice_ax.legend()



# integrate over solid angle, using heaviside weighing and different angles of acceptance
# cutoff angle
fig, axs = plt.subplots(1,2, layout='tight')

xstr = 'photon energies'
int_dict = interpolate_spectralrad_by_angle(cwv=10, x_vals=xstr, number_of_angles=100)
Ephs = int_dict[f'{xstr} (columns)']
dat_3D_int = int_dict['int rad']
angle_array = int_dict['angles (rows)']

Egs = np.arange(0.062, 0.20, 0.002)

cutoff_angles = np.arange(10,95,20)
for cutoff_angle in cutoff_angles:
    # downwelling photon flux (integrate over angles that are accepted)
    hvs = np.heaviside(cutoff_angle-angle_array, 1)
    rad_with_heaviside = np.transpose(np.transpose(dat_3D_int)*hvs)
    spectral_power_density = 2*np.trapz(rad_with_heaviside, np.radians(angle_array), axis=0)

    axs[0].plot(Ephs, spectral_power_density, label=cutoff_angle)

    phot_flux = spectral_power_density / (Ephs*q)  # [W.m-2/eV] / [eV][J/eV] --> [s-1.m-2/eV]
    dw_dict = {'photon energies':Ephs, 'downwelling photon flux':phot_flux}
    Ndot_vs_Eg = np.vectorize(Ndot_downwellheaviside)(Egs, dw_dict)
    axs[1].plot(Egs, Ndot_vs_Eg)

# diffusivity approx, for comparison
dw_diffusitivtyapprox = retrieve_downwelling_in_particleflux(cwv=10)
Ephs_diff = dw_diffusitivtyapprox['photon energies']
SPDF_diff = dw_diffusitivtyapprox['downwelling photon flux']*Ephs*q
axs[0].plot(Ephs_diff, SPDF_diff, ':k', label='Diffusivity approx (53$^\circ \\times \pi$)')
axs[1].plot(Egs, np.vectorize(Ndot_downwellheaviside)(Egs, dw_diffusitivtyapprox), ':k')

axs[0].legend()
axs[0].set_ylabel('Spectral Power Density Flux [W.m$^{-2}$/eV]')
axs[0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')

axs[1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
axs[1].set_xlabel('Bandgap, E$_g$ [eV]')



plt.show()
