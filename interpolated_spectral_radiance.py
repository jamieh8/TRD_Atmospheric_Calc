import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from TRD_Atmospheric_Functions import *
from matplotlib.colors import Normalize, LogNorm

def Eph_to_sx(x):
    return convert_from(x, units_in = 'photon energy [eV]', units_out = 'wavenumber [cm-1]')

def sx_to_Ephs(x):
    return convert_from(x, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')


def interpolation_method(x_predict, x_known, y_known, Eph):
    '''

    :param x_predict: Array of angles, in degrees, at which to interpolate data
    :param x_known: Array or list of angles at which data is known
    :param y_known: Data corresponding to known angles
    :return: Array of y values estimated for each of the x / angles requested
    '''

    # Transform x values for fitting/interpolation
    # p frac (p/p0 = 1/cos(theta))
    # x_predict = x_pred_split[0]
    new_x_predict = 1 / np.cos(np.radians(x_predict))
    new_x_known = 1 / np.cos(np.radians(x_known))

    # Interpolation + Extrapolation with scipy
    scp_int_ext = interpolate.interp1d(new_x_known, y_known, bounds_error=False, fill_value='extrapolate')
    y_predict = scp_int_ext(new_x_predict)

    return y_predict



# atm_data = atmospheric_dataset(cwv=24)
# cwv_str, Tskin = 'low', 301.56
cwv_str, Tskin = 'mid', 306.43
# cwv_str, Tskin = 'high', 299.86
atm_data = atmospheric_dataset_new(cwv=cwv_str, location='telfer', Tskin=Tskin)
emitter = planck_law_body(T=Tskin, Ephs=atm_data.photon_energies)

known_angles = atm_data.zenith_angles
col_heads = [f'downwelling_{theta}' for theta in known_angles]

Ephs = atm_data.photon_energies
wls = atm_data.wavelengths

angle_array = np.arange(0,91,1)


# make heatmap
include_slices = True
include_interpolation = True
plot_Lcos = True

plot_difference = False


ncols = 1
if include_slices:
    ncols += 1
if include_interpolation:
    ncols += 1


fig, axs = plt.subplots(1,ncols, layout='tight')

# dat_str = r'Spectral Radiance, L$_e$ [$\mathrm{W.m^{-2}.sr^{-1}/um}$]'
# xl_str = 'Wavelength [um]'
# xs = wls
# yvals = 'W.m-2'
# xvals = 'um'

# dat_str = r'Directional, Spectral PDF [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
dat_str = r'L$_\mathrm{ph}$ cos$\,\theta$ [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
# dat_str = r'$\Delta$ (L$_\mathrm{ph}$ cos$\,\theta$) [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
xl_str = 'Photon Energy, E$_\mathrm{ph}$ [eV]'
xs = Ephs
yvals = 's-1.m-2'
xvals = 'eV'


dat_2D = []
for col_head in col_heads:
    dat_2D += [atm_data.retrieve_spectral_array(yvals, xvals, col_head)]
dat_2D = np.array(dat_2D)  # in units [yvals.sr-1.xvals-1]

if plot_Lcos:
    dat_2D_int = atm_data.interpolate_cosDSPDF(angle_array)
else:
    dat_2D_int = atm_data.interpolate_by_angle(dat_2D, angle_array)

if plot_difference:
    # calculate difference between downwelling and upwelling
    # at each Eph & theta, calculate Lph costheta
    Lph_array = np.array([])
    for Eph in Ephs:
        Lph_array = np.append(Lph_array, emitter.angle_spectral_photon_flux(Eph=Eph, mu=-0.003))
    Lph_costheta = []
    for angle in angle_array:
        Lph_costheta += [Lph_array * np.cos(np.radians(angle))]
    Lph_costheta = np.array(Lph_costheta)

    dat_2D_int = Lph_costheta - dat_2D_int

    norm_0mid = Normalize(vmin=-1 * 1e22, vmax=1 * 1e22)
    style_args = {'cmap':'bwr', 'norm':norm_0mid, 'shading':'nearest'}
else:
    style_args = {'cmap': 'magma', 'norm':LogNorm(vmin=1e20, vmax=4e23), 'shading': 'nearest'}


if include_slices or include_interpolation:
    h_ax = axs[0]
else:
    h_ax = axs


hmap = h_ax.pcolor(xs, angle_array, dat_2D_int, **style_args)  # heatmap, 'inferno'
cbar = plt.colorbar(hmap)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(dat_str)
h_ax.set_ylabel('Zenith Angle, $\\theta$ [$\circ$]')
h_ax.set_xlabel(xl_str)


# secax = h_ax.secondary_xaxis('top', functions=(Eph_to_sx, sx_to_Ephs))
# secax.set_xlabel('Wavenumber [cm$^{-1}$]')
add_wl_ticks(h_ax)

h_ax.set_ylim([0,90])

# take slices of heatmap, plot separately as line plot
cmap = plt.get_cmap('tab10')
if include_slices:
    slice_ax = axs[1]

    # Wavelengths
    if xvals == 'um':
        slice_xs = [6,9.5,13,18]
        s_idxs = len(xs) - np.searchsorted(np.flip(xs), slice_xs)  # indexes of slice wavelengths
    else:
        slice_xs = [0.07, 0.08, 0.09, 0.1]
        s_idxs = np.searchsorted(Ephs, slice_xs)

    for i in range(len(s_idxs)):
        si = s_idxs[i]
        swav = slice_xs[i]
        Eph = atm_data.photon_energies[si]
        col = cmap(i)
        # plot spectral rad vs angle for a particular photon energy / wavelength
        spectral_rad_slice = dat_2D_int[:,si]
        slice_ax.plot(angle_array, spectral_rad_slice, '--', label=f'{swav} {xvals}', c=col)

        not_int_slice = dat_2D[:,si]  # not interpolated data at particular wavelength
        for angle_i, L_i in zip(known_angles, not_int_slice):
            # plot point for uninterpolated data
            y = L_i
            if plot_Lcos:
                y *= np.cos(np.radians(angle_i))
            slice_ax.plot(angle_i, y, 'o', c=col)

            if include_interpolation:
                axs[2].plot(1/np.cos(np.radians(angle_i)), L_i, 'o', c=col)  # plotted for each angle, each photon energy

        # angles_array_forlinint = np.arange(0,87,1)
        # new_x_interp = 1 / np.cos(np.radians(angles_array_forlinint))
        # new_x_known = 1 / np.cos(np.radians(known_angles))
        # interpolated_vals_forwav = interpolation_method(x_predict=new_x_interp, x_known=new_x_known, y_known=not_int_slice, Eph=Eph)
        # axs[2].plot(new_x_interp, interpolated_vals_forwav, '--', c=col)

        new_x = 1/np.cos(np.radians(angle_array[0:-3]))
        axs[2].plot(new_x, spectral_rad_slice[0:-3]*new_x, '--', c=col)
        h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--', c=col)

    slice_ax.set_xlabel('Zenith Angle, $\\theta$ [$\circ$]')
    slice_ax.set_ylabel(dat_str)
    slice_ax.legend()

if include_interpolation:
    axs[2].set_xlabel('Path Length Fraction, $\\frac{p(\\theta)}{p_0} = \\frac{1}{\cos\\theta}$')
    axs[2].set_ylabel(dat_str)
    axs[2].set_xlim([0.4,15])
    axs[2].set_ylabel(r'L$_\mathrm{ph}$ [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]')



plt.show()
