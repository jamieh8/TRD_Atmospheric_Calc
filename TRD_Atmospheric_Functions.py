import numpy as np
import os
import xarray as xr
from solcore.constants import kb, c, h, q
from scipy.optimize import minimize_scalar


def convert_from(array_in, units_in, units_out, corresponding_xs=None):
    if units_in == 'wavenumber [cm-1]':

        if units_out == 'wavelength [um]':
            return 1e4 / array_in
            # [1/cm-1]=[cm] ---> *1e4 --> [um]

        elif units_out == 'photon energy [eV]':
            phot_energies_J = h * c * (1e2*array_in)  # E = hf = hc v, [m-1]*[m.s-1][J.s] --> [J]
            return phot_energies_J / q  # [J] --> [eV]


    elif units_in == 'per wavenumber [/cm-1]':

        if units_out =='per wavelength [/um]':
            return array_in * corresponding_xs **2 * 1e-4
            # corresponding_xs = wavenumber [cm-1]
            # [Y/cm-1]*[cm-1]^2 = [Y/cm] --> *1e-4 --> [Y/um]

        elif units_out == 'per photon energy [/eV]':
            return 1e-2 * array_in / (h*c/q)
            # [cm-1/m-1][Y/cm-1]/[J.s][m.s-1][eV/J] --> [Y/m-1]/[eV.m] --> [Y/eV]


class atmospheric_dataset:
    def load_file_as_xarray(self, cwv):
        filename = os.path.join('simulations_telfer', f'telfer_australia_cwv{cwv}.txt')
        data = np.loadtxt(filename)

        # make an xarray to add headings etc.
        column_labels = [
            'wavenumber',  # cm-1, centre point of bin. 0.5 cm-1 steps
            'downwelling_0', 'downwelling_53', 'downwelling_70',  # RADIANCE - units of W cm-2 (cm-1)-1 sr-1
            'downwelling_flux', 'upwelling_flux', 'net_flux'  # FLUX - units of W cm-2 (cm-1)-1
        ]
        ds = xr.DataArray(data[:, 1:], coords={'wavenumber': data[:, 0],
                                               'column': column_labels[1:]}, dims=['wavenumber', 'column'])
        return ds

    def __init__(self, cwv):
        '''
        Loads in file as an xarray, with all original data values (i.e. no unit conversion), on initialization.
        Calculate wavenumbers, wavelengths, and photon_energies. Save resulting 1D arrays as class properties.
        :param cwv: column water vapour, should be 10, 24, or 54.
        '''
        self.org_xarray = self.load_file_as_xarray(cwv)

        self.wavenumbers = np.array(self.org_xarray['wavenumber'])
        self.wavelengths = convert_from(self.wavenumbers, units_in='wavenumber [cm-1]', units_out='wavelength [um]')
        self.photon_energies = convert_from(self.wavenumbers, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')

    def retrieve_spectral_array(self, yvals = 'W.m-2', xvals = 'eV', col_name = 'downwelling_flux'):
        '''
        :param yvals: Can be 'W.m-2', 's-1.m-2', or 'W.cm-2'. If string is not recognized, 'W.cm-2' is returned.
        :param xvals: Can be 'cm-1', 'eV', 'um'. If string is not recognized, 'cm-1' is returned.
        :param col_name: Can be 'downwelling_x' (with x = 0, 53, 70), 'downwelling_flux', 'upwelling_flux', 'net_flux'.
        :return: 1D array of spectral data specified by col_name, in units specified by yvals and xvals strings.
        '''
        new_array = self.org_xarray.sel(column=col_name)  # original values, in [W.cm-2.(sr-1)/cm-1]
        # (sr-1) depends on whether flux (no sr-1) or radiance (0,53,70) (with sr-1) is selected

        if yvals == 'W.m-2':
            # convert cm-2 to m-2
            new_array = new_array * 1e4  # [W.m-2.(sr-1)/cm-1]
        elif yvals == 's-1.m-2':
            new_array = new_array * 1e4 / (self.photon_energies * q)  #[W.m-2.(sr-1)/cm-1] / [eV][J/eV] --> [s-1.m-2.(sr-1)/cm-1]

        if xvals !='cm-1':
            if xvals == 'eV':
                units_out_str = 'per photon energy [/eV]'
            elif xvals == 'um':
                units_out_str = 'per wavelength [/um]'
            new_array = convert_from(new_array,
                                     units_in='per wavenumber [/cm-1]', units_out=units_out_str,
                                     corresponding_xs=self.wavenumbers)

        return new_array

    def interpolate_by_angle(self, dat_2D, angle_array):
        # dat_2D in [W.m-2/sr.eV] or [/sr.um]
        dat_2D_int = []
        for xi in range(len(dat_2D[0])):
            # for each wavelength, interpolate values for angles specified
            dat_2D_int += [np.interp(x=angle_array, xp=[0, 53, 70], fp=dat_2D[:, xi])]

        dat_2D_int = np.array(dat_2D_int)  # rows correspond to x values (photon energies or wavelengths)
        dat_2D_int = dat_2D_int.transpose()  # rows correspond to angles

        return dat_2D_int

    def interpolate_angle_spectral_data(self, angle_array, yvals = 'W.m-2', xvals='eV'):
        '''
        :param angle_array: 1D array of angles, in degrees, at which to interpolate the spectral data.
        :param yvals: 'W.m-2' or 's-1.m-2'
        :param xvals: 'cm-1', 'um', or 'eV'
        :return: 2D numpy array, with values in [yvals.sr-1.xvals-1]. Rows correspond to angles in angle_array. Columns corresponds to wavelength/photon energy/...
        '''
        # build 2D array with angles 0, 53, 70
        dat_2D = []
        for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
            dat_2D += [self.retrieve_spectral_array(yvals, xvals, col_head)]
        dat_2D = np.array(dat_2D)  # in units [yvals.sr-1.xvals-1]
        return self.interpolate_by_angle(dat_2D, angle_array)

    def retrieve_interpolated_angle_spectral_photflux(self, angle_array):
        '''
        Check if interpolated photon flux array exists (with right angle_array). If yes, return it. If no, run interpolation then return it.
        :param angle_array: 1D array of angles, in degrees, at which to interpolate the spectral data.
        :return: 2D numpy array, with values in [s-1.m-2.sr-1.eV-1]. Rows correspond to angles in angle_array. Columns corresponds to photon energies.
        '''

        if hasattr(self, 'int_ang_spec_photflux') and (self, 'angle_array') and np.array_equal(self.angle_array, angle_array):
            return self.int_ang_spec_photflux
        else:
            self.int_ang_spec_photflux = self.interpolate_angle_spectral_data(angle_array, yvals='s-1.m-2', xvals='eV')
            self.angle_array = angle_array
            return self.int_ang_spec_photflux

    def spectral_data_with_cutoffangle(self, angle_array, cutoff_angle):
        '''
        :param angle_array: 1D array of angles at which to interpolate spectral data, in degrees (0 to 90).
        :param cutoff_angle: between 0 and 90, in degrees. Angles smaller than cutoff will be "absorbed".
        :return: 1D array, spectral photon density flux [s-1.m-2/eV]
        '''
        photflux_int_2D = self.retrieve_interpolated_angle_spectral_photflux(angle_array)
        hvs = np.heaviside(cutoff_angle - angle_array, 1)

        #  method 1 - integrate over solid angles
        rad_with_heaviside = np.transpose(np.transpose(photflux_int_2D) * hvs)  # using heaviside to select "accepted" angles (*1), or "rejected" (*0)
        solid_angles = 2*np.pi*(1-np.cos(np.radians(angle_array)))
        spectral_photon_flux = np.trapz(rad_with_heaviside, solid_angles, axis=0)

        # method 2 - integrate over zenith angles
        # sin_zenith = np.sin(np.radians(angle_array))
        # rad_hv_sinz = 2*np.pi*np.transpose(np.transpose(photflux_int_2D) * hvs * sin_zenith)
        # spectral_photon_flux = np.trapz(rad_hv_sinz, np.radians(angle_array), axis=0)

        return spectral_photon_flux*0.5

    def retrieve_Ndot_heaviside(self, Eg, cutoff_angle):
        angle_array = np.linspace(0,90,100)
        Ephs = self.photon_energies
        spec_phot_flux = self.spectral_data_with_cutoffangle(angle_array, cutoff_angle)
        spec_pf_heavisided = spec_phot_flux * np.heaviside(Ephs-Eg, 0.5)
        return np.trapz(spec_pf_heavisided, Ephs)



class planck_law_body:
    def __init__(self, T=300):
        self.T = T
        self.kT_eV = kb * T / q  # [eV]

    def angle_spectral_photon_flux(self, Eph, mu):
        # [s-1.m-2 / sr.eV]
        return (2 / (c**2 * (h/q)**3)) * Eph**2 / (np.exp((Eph - mu) / self.kT_eV) - 1)

    def spectral_photon_flux(self, Eph, mu, cutoff_angle):
        # [s-1.m-2 / eV]
        sa = 2*np.pi*(1-np.cos(np.radians(cutoff_angle)))
        return 0.5 * self.angle_spectral_photon_flux(Eph, mu) * sa

    def retrieve_Ndot_heaviside(self, Ephs, Eg, mu, cutoff_angle):
        '''

        :param Ephs: Array of photon energies to use, in [eV]
        :param Eg: Bandgap, in [eV]
        :param mu: Fermi level split, in [eV]
        :return: Photon density flux [s-1.m-2]
        '''
        hvs_weight = np.heaviside(Ephs-Eg, 0.5)
        phot_flux = self.spectral_photon_flux(Ephs, mu, cutoff_angle)
        phot_flux_heavisided = hvs_weight*phot_flux

        return np.trapz(phot_flux_heavisided, Ephs)


class TRD_in_atmosphere:
    def __init__(self, TRD_body, atm_dataset, Ephs_TRD):
        self.TRD_body = TRD_body
        self.atm_data = atm_dataset
        self.Ephs = Ephs_TRD

    def power_density(self, mu, Eg, cutoff_angle):
        Ndot_out = self.TRD_body.retrieve_Ndot_heaviside(self.Ephs, Eg, mu, cutoff_angle)
        Ndot_in = self.atm_data.retrieve_Ndot_heaviside(Eg,cutoff_angle)
        J = q * (Ndot_out-Ndot_in)
        return J*mu

    def optimize_mu(self, Eg, cutoff_angle):
        opt_mu_dwh = minimize_scalar(self.power_density, bounds=[-Eg, 0],
                                     args=(Eg, cutoff_angle))
        return {'max power':opt_mu_dwh.fun, 'Vmpp':opt_mu_dwh.x}


def planck_dist(Eph, mu, kT):
    '''
    :param Eph: Particular photon energy, in eV
    :param mu: Fermi level splitting, in eV
    :param kT: in eV
    :return: Photon flux for the given eV, in [s-1.m-2.eV-1], calculated using generalized Planck's law (blackbody with temp T if mu=0).
    '''
    return ((2*np.pi)/(c**2*(h/q)**3))* Eph**2/(np.exp((Eph-mu)/kT)-1)

def spec_pflux_planckheaviside(Eph, Eg, mu, kT):
    '''
    Photon density flux calculated for using Planck's law for a semiconductor with Eg, at temp T and Fermi level split mu.
    Uses Heaviside weighing (i.e. 100% emission above Eg, 0 below Eg)
    :param Eph: Particular photon energy (or array) [eV]
    :param Eg: Bandgap [eV]
    :param mu: Fermi level splitting, [eV]
    :param kT: k*temperature of emitting body [eV]
    :return: Photon density flux [s-1.m-2] for the given eV. If Eph is an array, returns array of [s-1.m-2/eV] values.
    '''
    hvs_weight = np.heaviside(Eph-Eg, 0.5)
    pd_ys = planck_dist(Eph, mu, kT)
    return pd_ys*hvs_weight

def Ndot_boltzmann(Eg, T, mu): # particle flux density
    # accurate for large band gaps / large negative bias
    kT = kb*T/q # convert to eV to match units of Eg
    N = ((2*np.pi)/(c**2*h**3))*np.exp((mu-Eg)/kT)*kT*(Eg**2 + 2*Eg*kT + 2*kT**2)
    #eq 13 from Pusch et al 2019
    return N*q**3

def Ndot_planckheaviside(Ephs, Eg, mu, kT):
    spec_pflux = spec_pflux_planckheaviside(Ephs, Eg, mu, kT)
    pflux = np.trapz(spec_pflux, Ephs)
    return pflux

def Ndot_downwellheaviside(Eg, Ephs, downwell_array):
    get_heavisided = np.heaviside(Ephs - Eg, 0.5) * downwell_array
    pflux = np.trapz(get_heavisided, Ephs)
    return pflux

def neg_powerdensity_downwellheaviside(mu, Eg, Ephs_p, kT_converter, Ephs_atm, phot_flux_atm):
    # used to optimize over mu, with scipy
    Ndot_out = Ndot_planckheaviside(Ephs_p, Eg, mu, kT_converter)
    Ndot_in = Ndot_downwellheaviside(Eg, Ephs_atm, phot_flux_atm)
    J = q*(Ndot_out-Ndot_in)
    return mu*J

def neg_powerdensity_plancks(mu, Eg, Ephs, kT_convert, kT_env):
    # used to optimize over mu, with scipy
    Ndot_out = Ndot_planckheaviside(Ephs, Eg, mu, kT_convert)
    Ndot_in = Ndot_planckheaviside(Ephs, Eg, 0, kT_env)
    J = q * (Ndot_out - Ndot_in)
    return mu * J


