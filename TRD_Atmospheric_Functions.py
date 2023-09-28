import numpy as np
import os
import copy
import xarray as xr
import pygmo as pg
from solcore.constants import kb, c, h, q
from scipy.constants import sigma
from scipy.optimize import minimize_scalar
from scipy import interpolate, integrate


def convert_from(array_in, units_in, units_out, corresponding_xs=None):
    if units_in == 'wavenumber [cm-1]':

        if units_out == 'wavelength [um]':
            return 1e4 / array_in
            # [1/cm-1]=[cm] ---> *1e4 --> [um]

        elif units_out == 'photon energy [eV]':
            phot_energies_J = h * c * (1e2*array_in)  # E = hf = hc v, [m-1]*[m.s-1][J.s] --> [J]
            return phot_energies_J / q  # [J] --> [eV]

    if units_in == 'photon energy [eV]':
        if units_out == 'wavenumber [cm-1]':
            return array_in * q / (h*c*1e2)

        elif units_out == 'wavelength [um]':
            return 1e6 * (h/q) * c / array_in

    if units_in =='wavelength [um]':
        if units_out == 'photon energy [eV]':
            return 1e6 * (h/q) * c / array_in

    elif units_in == 'per wavenumber [/cm-1]':

        if units_out =='per wavelength [/um]':
            return array_in * corresponding_xs **2 * 1e-4
            # corresponding_xs = wavenumber [cm-1]
            # [Y/cm-1]*[cm-1]^2 = [Y/cm] --> *1e-4 --> [Y/um]

        elif units_out == 'per photon energy [/eV]':
            return 1e-2 * array_in / (h*c/q)
            # [cm-1/m-1][Y/cm-1]/[J.s][m.s-1][eV/J] --> [Y/m-1]/[eV.m] --> [Y/eV]


def Ndot_boltzmann(Eg, T, mu): # particle flux density
    # accurate for large band gaps / large negative bias
    kT = kb*T/q # convert to eV to match units of Eg
    N = ((2*np.pi)/(c**2*h**3))*np.exp((mu-Eg)/kT)*kT*(Eg**2 + 2*Eg*kT + 2*kT**2)
    #eq 13 from Pusch et al 2019
    return N*q**3

def Ndot_downwellheaviside(Eg, Ephs, downwell_array):
    get_heavisided = np.heaviside(Ephs - Eg, 0.5) * downwell_array
    pflux = np.trapz(get_heavisided, Ephs)
    return pflux


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
        :param cwv: column water vapour, should be 10, 24, or 54. Used to retrieve the correct file.
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
        elif yvals =='W.cm-2':
            pass
        else:
            print('yvals string not recognized. Returning spectral array in default W.cm-2')

        if xvals !='cm-1':
            if xvals == 'eV':
                units_out_str = 'per photon energy [/eV]'
            elif xvals == 'um':
                units_out_str = 'per wavelength [/um]'
            else:
                print('xvals string not recognized. Returning spectral array in default /eV')
                units_out_str = 'per photon energy [/eV]'
            new_array = convert_from(new_array,
                                     units_in='per wavenumber [/cm-1]', units_out=units_out_str,
                                     corresponding_xs=self.wavenumbers)

        return new_array

    def effective_skytemp(self, T_fill):
        downwelling_irrad_diff = self.retrieve_spectral_array(yvals='W.m-2', xvals='eV', col_name='downwelling_flux')
        Ephs = self.photon_energies
        downwelling_Wm2 = np.trapz(y=downwelling_irrad_diff, x=Ephs)

        # Ephs in dataset are limited. use blackbody with T_fill to "fill in" lower Ephs
        Ephs_fill = np.arange(0.0001,Ephs[0],Ephs[1]-Ephs[0])
        BB_fill = planck_law_body(T=T_fill, Ephs=Ephs_fill)
        BBfill_irrad = BB_fill.spectral_irradiance(Ephs_fill, mu=0, cutoff_angle=90)
        downwelling_Wm2 += np.trapz(y=BBfill_irrad, x=Ephs_fill)

        effectiveT = (downwelling_Wm2/sigma)**(0.25)   # ( [W.m-2]/[W.m-2.K-4] ) ^ 1/4 = [K]
        # sigma : stefan-boltzmann constant

        return effectiveT

    def interpolate_by_angle(self, dat_2D, angle_array):
        '''
        :param dat_2D:
        :param angle_array: Array of angles (in degrees) at which to interpolate
        :return: Interpolated 2D array, in same using as dat_2D. Rows correspond to angles, columns to same spectral units as dat_2D.
        '''
        # dat_2D in [W.m-2/sr.eV] or [/sr.um]
        dat_2D_int = []
        for xi in range(len(dat_2D[0])):
            # for each wavelength/photon energy, interpolate values for angles specified
            new_x_predict = 1 / np.cos(np.radians(angle_array))
            new_x_known = 1 / np.cos(np.radians([0,53,70]))
            y_known = dat_2D[:, xi]
            scp_int_ext = interpolate.interp1d(new_x_known, y_known, bounds_error=False, fill_value='extrapolate')
            y_predict = scp_int_ext(new_x_predict)

            dat_2D_int += [y_predict]

        dat_2D_int = np.array(dat_2D_int)  # rows correspond to x values (photon energies or wavelengths)
        dat_2D_int = dat_2D_int.transpose()  # rows correspond to angles

        return dat_2D_int

    def interpolate_angle_spectral_data(self, angle_array, yvals = 'W.m-2', xvals='eV'):
        '''
        Retrieves 2D array with directional spectral radiance (expressed in y units requested, plus 'sr-1' for directionality).
        Passes this 2D array to interpolate_by_angle, to interpolate according to the angle array requested.
        :param angle_array: 1D array of angles, in degrees, at which to interpolate the spectral data.
        :param yvals: 'W.m-2' or 's-1.m-2'
        :param xvals: 'cm-1', 'um', or 'eV'
        :return: 2D numpy array, with values in [yvals.sr-1.xvals-1]. Rows correspond to angles in angle_array. Columns corresponds to spectral values.
        '''
        # build 2D array with angles 0, 53, 70
        dat_2D = []
        for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
            dat_2D += [self.retrieve_spectral_array(yvals, xvals, col_head)]
        dat_2D = np.array(dat_2D)  # in units [yvals.sr-1.xvals-1]
        return self.interpolate_by_angle(dat_2D, angle_array)

    def interpolate_cosDSPDF(self, angle_array):
        '''
        :param angle_array: Angle at which to interpolate Lph*costheta
        :return: L_ph*cos(theta), calculated by interpolating (and extrapolating, where needed) L_2D at the angles requested
        '''

        dat_2D_interpolated = []

        angles_rad = np.radians(angle_array)
        angles_known_rad = np.radians([0,53,70])
        L_2D = self.Lph_2D  # 2D array of L_ph vals, containing directional spectral photon flux [s-1.m-2.sr-1/eV] at 0, 53, 70 degrees

        for theta in angles_rad:
            # check which pair of angles to use
            # usually (theta1 < theta < theta2) for interpolation. theta2 < theta for extrapolation.
            if theta < angles_known_rad[1]:
                it1, it2 = 0, 1
            else:
                it1, it2 = 1, 2
            theta1, theta2 = angles_known_rad[it1], angles_known_rad[it2]

            angle_multiplier = (1-np.cos(theta)/np.cos(theta1)) / (1/np.cos(theta2)-1/np.cos(theta1))

            L_costheta_int = []
            for xi in range(len(L_2D[0])):
                # for each wavelength / photon energy
                L_t1 = L_2D[it1][xi]
                L_t2 = L_2D[it2][xi]

                L_costheta = (L_t2-L_t1)*angle_multiplier + L_t1*np.cos(theta)
                L_costheta_int += [L_costheta]

            dat_2D_interpolated += [L_costheta_int]

        return np.array(dat_2D_interpolated)

    def retrieve_interpolated_cosDSPDF(self, angle_array):
        '''
        Check if interpolated photon flux array exists (with right angle_array). If yes, return it. If no, run interpolation then return it.
        :param angle_array: 1D array of zenith angles, in degrees, at which to interpolate the spectral data.
        :return: 2D numpy array, with values in [s-1.m-2.sr-1.eV-1]. Rows correspond to the zenith angles in angle_array. Columns corresponds to photon energies.
        '''

        if hasattr(self, 'interpolated_cosDSPDF') and hasattr(self, 'angle_array') and np.array_equal(self.angle_array, angle_array):
            # if interpolated array already exists, and corresponds to angle_array requested, return existing data
            return self.interpolated_cosDSPDF
        else:
            # populate Lph_2D
            dat_2D = []
            for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
                dat_2D += [self.retrieve_spectral_array('s-1.m-2', 'eV', col_head)]
            dat_2D = np.array(dat_2D)  # in units [s-1.m-2.sr-1/eV]
            self.Lph_2D = dat_2D

            # interpolate across zenith angle and store results
            self.interpolated_cosDSPDF = self.interpolate_cosDSPDF(angle_array)
            self.angle_array = angle_array
            return self.interpolated_cosDSPDF

    def spectral_PDF_with_cutoffangle(self, angle_array, cutoff_angle):
        '''
        :param angle_array: 1D array of angles at which to interpolate spectral data, in degrees (0 to 90).
        :param cutoff_angle: zenith angle between 0 and 90, in degrees. Angles smaller than or equal to the cutoff will be 100% "absorbed".
        :return: 1D array, spectral photon density flux [s-1.m-2/eV]
        '''
        # zenith angle theta, azimuth angle phi
        # d(Omega) = sin(theta) * d(theta) * d(phi), element solid angle
        # cos(theta) factor for Lambertian emitter
        # Spectral_PDF = \int_0^cutoff Directional_Spectral_PDF cos(theta) d(Omega)
        # Spectral_PDF = \int_0^2pi d(phi) \int_0^theta_cutoff Directional_Spectral_PDF cos(theta) sin(theta) d(theta)
        # Spectral_PDF = 2pi \int_0^90 heaviside * Directional_Spectral_PDF cos(theta) sin(theta) d(theta)

        # Using pre-existing sampled array
        # DSPDF * cos(theta) is computed altogether such to avoid ill-defined DSPDF at theta near 90..
        DSPDF_cos = self.retrieve_interpolated_cosDSPDF(angle_array)
        hvs = np.heaviside(cutoff_angle - angle_array, 1)  # using heaviside to select "accepted" angles (*1), or "rejected" (*0)
        sin_zenith = np.sin(np.radians(angle_array))  # zenith angle -> theta
        rad_hv_sinz = np.transpose(np.transpose(DSPDF_cos) * hvs * sin_zenith)  # [s-1.m-2.sr-1/eV]*[rad]

        integral_over_theta = np.trapz(rad_hv_sinz, np.radians(angle_array), axis=0)

        spectral_photon_flux = 2 * np.pi *  integral_over_theta  # [s-1.m-2.rad-1/eV]*[rad]
        return spectral_photon_flux  # [s-1.m-2/eV]


    def retrieve_Ndot_heaviside(self, Eg, cutoff_angle):
        if cutoff_angle == None:
            # if cutoff_angle is None, use the diffusivity approximation
            spec_phot_flux = self.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV', col_name='downwelling_flux')
        else:
            # if cutoff angle is given, perform integral over interpolated angles
            angle_array = np.arange(0,90.5,0.5)
            spec_phot_flux = self.spectral_PDF_with_cutoffangle(angle_array, cutoff_angle)

        Ephs = self.photon_energies
        spec_pf_heavisided = spec_phot_flux * np.heaviside(Ephs-Eg, 0.5)
        int_over_Eph = np.trapz(spec_pf_heavisided, Ephs)

        return int_over_Eph



class planck_law_body:
    def __init__(self, T=300, Ephs=np.arange(1e-6, 0.31, 0.0001)):
        self.T = T
        self.kT_eV = kb * T / q  # [eV]
        self.Ephs = Ephs

    def angle_spectral_photon_flux(self, Eph, mu):
        '''
        :param Eph: Eph in eV
        :param mu: mu in eV
        :return: spectral, directional photon flux - value, or array of values, in [s-1.m-2.sr-1/eV]
        '''
        # [s-1.m-2 / sr.eV]
        return (2 / (c**2 * (h/q)**3)) * Eph**2 / (np.exp((Eph - mu) / self.kT_eV) - 1)

    def spectral_photon_flux(self, Eph, mu=0, cutoff_angle = 90):
        '''
        :param Eph: Photon energy/energies, in [eV]
        :param mu: Fermi level splitting, in [eV]
        :param cutoff_angle: Angle between 0 and 90 [degrees], past which no emission occurs.
        :return: Spectral photon flux - value or array of values in [s-1.m-2/eV]
        '''
        if cutoff_angle == None:
            cutoff_angle = 90
        # [s-1.m-2.eV-1] = [s-1.m-2.sr-1.eV-1]*[sr] -- must integrate over solid angle
        # d(Omega) = sin(theta) * d(theta) * d(phi), element solid angle for int
        # angle_spectral_photon_flux is independent of theta, can be taken out.
        # Spectral_PDF = \int_0^cutoff Directional_Spectral_PDF cos(theta) d(Omega)
        # Spectral_PDF = \int_0^2pi d(phi) \int_0^theta_cutoff Directional_Spectral_PDF cos(theta) sin(theta) d(theta)
        # Spectral_PDF = 2pi * Directional_Spectral_PDF * \int_0^theta_cutoff cos(theta) sin(theta) d(theta)
        # with \int_0^theta_cutoff cos(theta) sin(theta) d(theta) = [sin^2(theta)/2]_0^theta_cutoff = sin^2(theta_cutoff)/2
        int_sinz_cosz = np.sin(np.radians(cutoff_angle))**2 / 2
        return 2*np.pi * self.angle_spectral_photon_flux(Eph, mu) * int_sinz_cosz

    def spectral_irradiance(self, Eph, mu=0, cutoff_angle=90):
        '''
        :param Eph:
        :param mu:
        :param cutoff_angle:
        :return: value or array of values in [W.m-2/eV]
        '''
        return self.spectral_photon_flux(Eph, mu, cutoff_angle) * q * Eph  # [s-1.m-2/eV]*[J/eV]*[eV]

    def retrieve_Ndot_heaviside(self, Eg, cutoff_angle=90, mu=0):
        '''
        :param Eg: Bandgap, in [eV]
        :param mu: Fermi level split, in [eV]
        :return: Photon density flux [s-1.m-2]
        '''
        Ephs = self.Ephs
        hvs_weight = np.heaviside(Ephs-Eg, 0.5)
        if cutoff_angle == None:
            cutoff_angle = 90
            # cutoff angle of None used for diffusivity approx. in planck body case, diff = 90 (same "cost")
        phot_flux = self.spectral_photon_flux(Ephs, mu, cutoff_angle)
        phot_flux_heavisided = hvs_weight*phot_flux

        return np.trapz(phot_flux_heavisided, Ephs)


class TRD_in_atmosphere:
    def __init__(self, TRD_body, atmosphere):
        self.TRD_body = TRD_body
        self.atm = atmosphere

    def current_density(self, mu, Eg, cutoff_angle, eta_ext=1, consider_nonrad = False):
        Ndot_out = self.TRD_body.retrieve_Ndot_heaviside(Eg = Eg, mu = mu, cutoff_angle = cutoff_angle)
        Ndot_in = self.atm.retrieve_Ndot_heaviside(Eg = Eg, cutoff_angle = cutoff_angle)

        if consider_nonrad:
            Ndot_out_mu0 = self.TRD_body.retrieve_Ndot_heaviside(Eg = Eg, mu = 0, cutoff_angle = cutoff_angle)
            J = q * (Ndot_out/eta_ext - Ndot_in + Ndot_out_mu0*(1-1/eta_ext))
        else:
            J = q * (Ndot_out - Ndot_in)

        return J

    def power_density(self, mu, Eg, cutoff_angle, eta_ext=1, consider_nonrad=False):
        J = self.current_density(mu, Eg, cutoff_angle, eta_ext, consider_nonrad)
        return J * mu

    def optimize_mu(self, Eg, cutoff_angle):
        opt_mu_dwh = minimize_scalar(self.power_density, bounds=[-Eg, 0],
                                     args=(Eg, cutoff_angle))
        return {'max power':opt_mu_dwh.fun, 'Vmpp':opt_mu_dwh.x}


class optimize_powerdensity:
    def __init__(self, trd_in_environment, args_to_opt, args_to_fix):
        self.trd_in_env = trd_in_environment
        self.args_to_fix = args_to_fix  # dict with fixed values to use
        self.args_to_opt = args_to_opt  # list of strings with arguments to optimize

    def fitness(self, x):
        # assemble arguments
        new_args = self.args_to_fix
        for ai, opt_args in enumerate(self.args_to_opt):
            if opt_args == 'mu_frac':
                new_args.update({'mu':x[ai]*new_args['Eg']})
                # ^ require 'mu_frac' to be populated last, so that 'Eg' is already defined
            else:
                new_args.update({opt_args:x[ai]})


        power_density = self.trd_in_env.power_density(**new_args)
        return [power_density]

    def get_bounds(self):
        # bounds for each argument, for reference
        if isinstance(self.trd_in_env.atm, atmospheric_dataset):
            Eg_min = 0.062  # limit min Eg to smallest photon energy in data
        else:
            Eg_min = 1e-4
        bound_ref = {'Eg':{'min':Eg_min, 'max':0.15},
                     'mu_frac':{'min':-1, 'max':0},
                     'cutoff_angle':{'min':10, 'max':90},
                     'eta_ext':{'min':0.001, 'max':1}}

        # retrieve relevant bounds, for the arguments being optimized
        mins, maxs = [], []
        for opt_args in self.args_to_opt:
            mins += [bound_ref[opt_args]['min']]
            maxs += [bound_ref[opt_args]['max']]

        return (mins, maxs)

    def get_nix(self):
        # integer dimension
        if 'cutoff_angle' in self.args_to_opt:
            return 1
        else:
            return 0


def get_best_pd(trd_in_environment, args_to_opt, args_to_fix, alg):
    '''
    :param trd_in_environment: Instance of TRD_in_atmosphere class, with method power_density defined
    :param args_to_opt: List of strings defining which arguments of power_density to optimize. Options are 'Eg', 'mu', 'cutoff_angle'. Order of strings determines order in which optimized x values are returned.
    :param args_to_fix: Dictionary, with keys that are strings specifying the arguments ('Eg', 'mu', 'cutoff_angle') to fix for the optimization. Corresponding values are values to fix.
    :param alg: pygmo algorithm to use for optimization
    :return: best_x (dictionary with keys corresponding to args_to_opt strings, values giving parameters corresponding to best_pd), best_pd ("smallest" power density found, W.m-2)
    '''


    args_to_opt_mod = copy.deepcopy(args_to_opt)
    if 'mu' in args_to_opt_mod:
        # if mu is to be optimized, replace ID string with 'mu_frac', and ensure it is last
        args_to_opt_mod[args_to_opt.index('mu')] = 'mu_frac'

    if 'cutoff_angle' in args_to_opt_mod:
        # if cutoff_angle is to be optimized, ensure it is at the end of the list (requirement for pygmo integer limit)
        args_to_opt_mod.remove('cutoff_angle')
        args_to_opt_mod += ['cutoff_angle']

    # define problem with UDP
    prob = pg.problem(optimize_powerdensity(trd_in_environment, args_to_opt_mod, args_to_fix))

    # initial population
    pop = pg.population(prob, size=15)
    algo = pg.algorithm(alg)
    # algo.set_verbosity(1)   # toggle for debugging

    pop = algo.evolve(pop)
    best_pd = pop.champion_f
    best_x = pop.champion_x

    x_opt_dict = {}
    for bxi, argo in zip(best_x, args_to_opt_mod):
        if argo == 'mu_frac':
            # get Eg
            try:
                Eg = args_to_fix['Eg']
            except:
                Eg = x_opt_dict['Eg']
            x_opt_dict.update({'mu':Eg*bxi})
        else:
            x_opt_dict.update({argo:bxi})

    return x_opt_dict, best_pd
