import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

def Eph_to_sx(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def sx_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')



atm_data = atmospheric_dataset(cwv=10)

Ephs = atm_data.photon_energies
wls = atm_data.wavelengths

angle_array = np.linspace(0,90,10)


# make heatmap
include_slices = True
wls_and_rad = True

if include_slices:
    ncols = 2
else:
    ncols = 1

fig, axs = plt.subplots(1,ncols, layout='tight')

if wls_and_rad:
    dat_3D_int = atm_data.interpolate_angle_spectral_data(angle_array, yvals='W.m-2', xvals='um')
    dat_str = r'Spectral Radiance [$\mathrm{W.m^{-2}.sr^{-1}/um}$]'
    xl_str = 'Wavelength [um]'
    xs = wls
else:
    dat_3D_int = atm_data.interpolate_angle_spectral_data(angle_array, yvals='s-1.m-2', xvals='eV')
    dat_str = r'Directional, Spectral PDF [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
    xl_str = 'Photon Energy, E$_\mathrm{ph}$ [eV]'
    xs = Ephs


if include_slices:
    h_ax = axs[0]
else:
    h_ax = axs
hmap = h_ax.pcolor(xs, angle_array, dat_3D_int, cmap='plasma', shading='nearest')  # heatmap
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(dat_str)
h_ax.set_ylabel('Zenith Angle, $\\theta$ [$\circ$]')
h_ax.set_xlabel(xl_str)

# secax = h_ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
# secax.set_xlabel('Wavenumber [cm$^{-1}$]')


# take slices of heatmap, plot separately as line plot
if include_slices:
    slice_ax = axs[1]

    # Wavelengths
    if wls_and_rad:
        slice_xs = [6,9.5,13,18]
        s_idxs = len(xs) - np.searchsorted(np.flip(xs), slice_xs)  # indexes of slice wavelengths
    else:
        slice_xs = [0.08, 0.12, 0.15, 0.2]
        s_idxs = np.searchsorted(Ephs, slice_xs)

    for si, swav in zip(s_idxs, slice_xs):
        # plot spectral rad vs angle for a particular photon energy / wavelength
        spectral_rad_slice = dat_3D_int[:,si]
        slice_ax.plot(angle_array, spectral_rad_slice, label=f'{swav}')

        h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--')

    slice_ax.set_xlabel('Zenith Angle, $\\theta$ [$\circ$]')
    slice_ax.set_ylabel(dat_str)
    slice_ax.legend()



# integrate over solid angle, using heaviside weighing and different angles of acceptance
# cutoff angle
# fig, axs = plt.subplots(1,2, layout='tight')
#
# Egs = np.arange(0.062, 0.20, 0.002)
#
# cutoff_angles = np.arange(10,95,20)
# for cutoff_angle in cutoff_angles:
#     # downwelling photon flux (integrate over angles that are accepted)
#     spectral_photon_flux = atm_data.spectral_data_with_cutoffangle(angle_array, cutoff_angle)
#     axs[0].plot(Ephs, spectral_photon_flux, label=f'cutoff {cutoff_angle}$^\circ$')
#
#     Ndot_vs_Eg = np.vectorize(atm_data.retrieve_Ndot_heaviside)(Egs, cutoff_angle)
#     axs[1].plot(Egs, Ndot_vs_Eg)
#
#
# # diffusivity approx, for comparison
# SPDF_diff = atm_data.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV', col_name='downwelling_flux')
# Ephs_diff = atm_data.photon_energies
# axs[0].plot(Ephs_diff, SPDF_diff, ':k', label='Diffusivity approx (53$^\circ \\times \pi$)')
# axs[1].plot(Egs, np.vectorize(Ndot_downwellheaviside, excluded=[1,2])(Egs, Ephs_diff, SPDF_diff), ':k')
#
#
# axs[0].legend()
# axs[0].set_ylabel('Spectral Photon Density Flux [s$^{-1}$.m$^{-2}$/eV]')
# axs[0].set_xlabel('Photon Energy, E$_{ph}$ [eV]')
#
# axs[1].set_ylabel('Photon Density Flux, $\dot{\mathrm{N}}$ [m$^{-2}$.s$^{-1}]$')
# axs[1].set_xlabel('Bandgap, E$_g$ [eV]')
#
# for ax in axs:
#     secax = ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
#     secax.set_xlabel('Wavenumber [cm$^{-1}$]')

plt.show()
