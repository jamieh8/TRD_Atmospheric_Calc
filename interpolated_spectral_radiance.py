import matplotlib.pyplot as plt
from scipy import interpolate
from TRD_Atmospheric_Functions import *

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

    # Chop out angles > 80 for separate handling
    # ind_80 = np.where(x_predict==80)[0]
    # x_pred_split = np.split(x_predict, ind_80)
    # x_predict_geq80 = x_pred_split[1]
    # bb = planck_law_body(T=300, Ephs=Eph)
    # SDphotflux = bb.angle_spectral_photon_flux(Eph, 0)
    # rad_per_eV = q * Eph * SDphotflux  # [s-1.m-2.sr-1/eV]*[eV]*[J/eV]
    # rad_per_um = rad_per_eV * Eph ** 2 / (h * c / q) * 1e-6  # [W.m-2.sr-1/eV] * [eV] / [um]
    #^ use this value for angles > 80

    # Transform x values for fitting/interpolation
    # p frac (p/p0 = 1/cos(theta))
    # x_predict = x_pred_split[0]
    new_x_predict = 1 / np.cos(np.radians(x_predict))
    new_x_known = 1 / np.cos(np.radians(x_known))

    # x = p-p0 = p0*(1/costheta - 1)
    # new_x_known = 1 / np.cos(np.radians(x_known)) - 1
    # new_x_predict = 1 / np.cos(np.radians(x_predict)) - 1

    # Basic interpolation with numpy
    # y_predict = np.interp(x=x_predict, xp=x_known, fp=y_known)

    # Linear Fit with numpy
    # fit_parms = np.polyfit(x=x_known, y=y_known, deg=1)
    # y_predict = fit_parms[0]*x_predict + fit_parms[1]

    # Interpolation + Extrapolation with scipy
    scp_int_ext = interpolate.interp1d(new_x_known, y_known, bounds_error=False, fill_value='extrapolate')
    y_predict = scp_int_ext(new_x_predict)

    return y_predict



atm_data = atmospheric_dataset(cwv=10)

Ephs = atm_data.photon_energies
wls = atm_data.wavelengths

angle_array = np.arange(0,90,1)


# make heatmap
include_slices = True
include_interpolation = True
plot_Lcos = True

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

dat_str = r'Directional, Spectral PDF [$\mathrm{s^{-1}.m^{-2}.sr^{-1}/eV}$]'
xl_str = 'Photon Energy, E$_\mathrm{ph}$ [eV]'
xs = Ephs
yvals = 's-1.m-2'
xvals = 'eV'


dat_2D = []
for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
    dat_2D += [atm_data.retrieve_spectral_array(yvals, xvals, col_head)]
dat_2D = np.array(dat_2D)  # in units [yvals.sr-1.xvals-1]

if plot_Lcos:
    dat_2D_int = atm_data.interpolated_cosDSPDF(angle_array)
else:
    dat_2D_int = atm_data.interpolate_by_angle(dat_2D, angle_array)


if include_slices or include_interpolation:
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
    if xvals == 'um':
        slice_xs = [6,9.5,13,18]
        s_idxs = len(xs) - np.searchsorted(np.flip(xs), slice_xs)  # indexes of slice wavelengths
    else:
        slice_xs = [0.07, 0.08, 0.09, 0.1]
        s_idxs = np.searchsorted(Ephs, slice_xs)

    for i in range(len(s_idxs)):
        si = s_idxs[i]
        swav = slice_xs[i]
        col = cmap(i)
        # plot spectral rad vs angle for a particular photon energy / wavelength
        spectral_rad_slice = dat_2D_int[:,si]
        slice_ax.plot(angle_array, spectral_rad_slice, '--', label=f'{swav} {xvals}', c=col)

        not_int_slice = dat_2D[:,si]  # not interpolated data at particular wavelength
        for angle_i, L_i in zip([0,53,70], not_int_slice):

            # plot point for uninterpolated data
            y = L_i
            if plot_Lcos:
                y *= np.cos(np.radians(angle_i))
            slice_ax.plot(angle_i, y, 'o', c=col)

            if include_interpolation:
                Eph = atm_data.photon_energies[si]

                # new_x_ref = 1/np.cos(np.radians([0,53,70]))
                # y_ref = not_int_slice * np.cos(np.radians([0,53,70]))
                axs[2].plot(1/np.cos(np.radians(angle_i)), y, 'o', c=col)

                new_x_interp = 1 / np.cos(np.radians(angle_array))
                interpolated_vals_forwav = interpolation_method(x_predict=angle_array, x_known=[0,53,70], y_known = not_int_slice, Eph= Eph)
                axs[2].plot(new_x_interp, interpolated_vals_forwav, '--', c=col)

        h_ax.plot(2*[swav], [angle_array[0],angle_array[-1]], '--', c=col)

    slice_ax.set_xlabel('Zenith Angle, $\\theta$ [$\circ$]')
    slice_ax.set_ylabel(dat_str)
    slice_ax.legend()

if include_interpolation:
    axs[2].set_xlabel('Path Length Fraction, $\\frac{p(\\theta)}{p_0} = \\frac{1}{\cos\\theta}$')
    axs[2].set_ylabel(dat_str)

plt.show()
