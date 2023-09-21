import matplotlib.pyplot as plt
from scipy import interpolate
from TRD_Atmospheric_Functions import *

def Eph_to_sx(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def sx_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')

def interpolation_method(x_predict, x_known, y_known):
    # y_predict = np.interp(x=x_predict, xp=x_known, fp=y_known)

    # fit_parms = np.polyfit(x=x_known, y=y_known, deg=1)
    # y_predict = fit_parms[0]*x_predict + fit_parms[1]

    scp_int_ext = interpolate.interp1d(x_known, y_known, bounds_error=False, fill_value='extrapolate')
    y_predict = scp_int_ext(x_predict)
    return y_predict

atm_data = atmospheric_dataset(cwv=10)

Ephs = atm_data.photon_energies
wls = atm_data.wavelengths

angle_array = np.arange(0,90,5)


# make heatmap
include_slices = True
include_interpolation = True
wls_and_rad = True

ncols = 1
if include_slices:
    ncols += 1
if include_interpolation:
    ncols += 1


fig, axs = plt.subplots(1,ncols, layout='tight')

# if wls_and_rad:
#     dat_2D_int = atm_data.interpolate_angle_spectral_data(angle_array, yvals='W.m-2', xvals='um')
#     dat_str = r'Spectral Radiance, L$_e$ [$\mathrm{W.m^{-2}.sr^{-1}/um}$]'
#     xl_str = 'Wavelength [um]'
#     xs = wls
# else:
#     dat_2D_int = atm_data.interpolate_angle_spectral_data(angle_array, yvals='s-1.m-2', xvals='eV')
#     dat_str = r'Directional, Spectral PDF [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
#     xl_str = 'Photon Energy, E$_\mathrm{ph}$ [eV]'
#     xs = Ephs

dat_str = r'Spectral Radiance, L$_e$ [$\mathrm{W.m^{-2}.sr^{-1}/um}$]'
xl_str = 'Wavelength [um]'
xs = wls
yvals = 'W.m-2'
xvals = 'um'
dat_2D = []
for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
    dat_2D += [atm_data.retrieve_spectral_array(yvals, xvals, col_head)]
dat_2D = np.array(dat_2D)  # in units [yvals.sr-1.xvals-1]

dat_2D_int = []
for xi in range(len(dat_2D[0])):
    # for each wavelength, interpolate values for angles specified
    Le_wav = dat_2D[:,xi] # L_e values for this wavelength
    pfrac_vals = 1/np.cos(np.radians(angle_array))
    pfrac_ref = 1/np.cos(np.radians([0,53,70]))
    interpolated_vals_forwav = interpolation_method(pfrac_vals, pfrac_ref, Le_wav)
    dat_2D_int += [interpolated_vals_forwav]
dat_2D_int = np.array(dat_2D_int)  # rows correspond to x values (photon energies or wavelengths)
dat_2D_int = dat_2D_int.transpose()


if include_slices:
    h_ax = axs[0]
else:
    h_ax = axs
hmap = h_ax.pcolor(xs, angle_array, dat_2D_int, cmap='inferno', shading='nearest')  # heatmap
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(dat_str)
h_ax.set_ylabel('Zenith Angle, $\\theta$ [$\circ$]')
h_ax.set_xlabel(xl_str)



# secax = h_ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
# secax.set_xlabel('Wavenumber [cm$^{-1}$]')


# take slices of heatmap, plot separately as line plot
cmap = plt.get_cmap('tab10')
if include_slices:
    slice_ax = axs[1]

    # Wavelengths
    if wls_and_rad:
        slice_xs = [6,9.5,13,18]
        s_idxs = len(xs) - np.searchsorted(np.flip(xs), slice_xs)  # indexes of slice wavelengths
    else:
        slice_xs = [0.08, 0.12, 0.15, 0.2]
        s_idxs = np.searchsorted(Ephs, slice_xs)

    for i in range(len(s_idxs)):
        si = s_idxs[i]
        swav = slice_xs[i]
        col = cmap(i)
        # plot spectral rad vs angle for a particular photon energy / wavelength
        spectral_rad_slice = dat_2D_int[:,si]
        slice_ax.plot(angle_array, spectral_rad_slice, label=f'{swav}', c=col)


        not_int_slice = dat_2D[:,si]  # not interpolated data at particular wavelength
        for angle_i, L_i in zip([0,53,70], not_int_slice):
            slice_ax.plot(angle_i, L_i, 'o', c=col)

            if include_interpolation:
                pfrac_ref = 1/np.cos(np.radians([0,53,70]))
                axs[2].plot(pfrac_ref, not_int_slice, 'o', c=col)
                pfrac_vals = 1 / np.cos(np.radians(angle_array))
                interpolated_vals_forwav = interpolation_method(pfrac_vals, pfrac_ref, not_int_slice)
                axs[2].plot(pfrac_vals, interpolated_vals_forwav, '--', c=col)

        h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--', c=col)

    slice_ax.set_xlabel('Zenith Angle, $\\theta$ [$\circ$]')
    slice_ax.set_ylabel(dat_str)
    slice_ax.legend()

if include_interpolation:
    axs[2].set_xlabel('Path($\\theta$)/Path($\\theta$=0), $\\frac{p}{a}$')
    axs[2].set_ylabel(dat_str)

plt.show()
