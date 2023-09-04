import numpy as np
import os
import xarray as xr
from solcore.constants import kb, c, h, q

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

def load_file_as_xarray(cwv = 10):#, convert_to_wavelength = True):
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

    # if convert_to_wavelength:
    #     # converting wavenumber to wavelength: just 1/wavenumber. But also need to adjust flux/radiance since these
    #     # are per cm-1 (Jacobian transformation)
    #     ds = 1e-4 * ds * (ds['wavenumber']) ** 2  # [W.cm-2/cm-1]*[cm-1]^2 = [W.cm-2/cm] --> *1e-4 --> [W.cm-2/um]
    #     ds = ds.rename({'wavenumber': 'wavelength'})
    #     ds.coords['wavelength'] = 1e4 / ds.coords['wavelength']  # [1/cm-1]=[cm] ---> *1e4 --> [um]

    return ds


def retrieve_downwellingrad_as_nparray(cwv=10, x_vals = 'photon energies'):
    '''

    :param cwv: column water vapour. 10, 24, or 54
    :return: dict with keys 'wavelengths' (1D array of wavelengths in um) and 'downwelling rad' (2D array in [W.m-2.sr-1/um]). Rows are degrees. Columns correspond to wavelengths.
    '''
    dat_xarray = load_file_as_xarray(cwv)
    wavnums = dat_xarray['wavenumber']  # [cm-1]

    dat_3D = []
    for col_head in ['downwelling_0', 'downwelling_53', 'downwelling_70']:
        dat_3D += [dat_xarray.sel(column=col_head)]
    dat_3D = np.array(dat_3D)  # [W.m-2.sr-1/um]

    if x_vals == 'wavelengths':
        x_array = convert_from(wavnums, units_in='wavenumber [cm-1]', units_out='wavelength [um]')
        dat_3D_per_um = convert_from(dat_3D, units_in='per wavenumber [/cm-1]', units_out='per wavelength [/um]',
                                     corresponding_xs=np.array(3*[wavnums]))   # [W.cm-2.sr-1/cm-1] --> [W.cm-2.sr-1/um]
        dat_3D_converted = 1e4 * dat_3D_per_um  # [W.m-2.sr-1/um]

    elif x_vals == 'photon energies':
        x_array = convert_from(wavnums, units_in='wavenumber [cm-1]', units_out='photon energy [eV]')
        dat_3D_per_um = convert_from(dat_3D, units_in='per wavenumber [/cm-1]', units_out='per photon energy [/eV]')   # [W.cm-2.sr-1/cm-1] --> [W.cm-2.sr-1/eV]
        dat_3D_converted = 1e4 * dat_3D_per_um  # [W.m-2.sr-1/um]

    return {x_vals:x_array, 'downwelling rad':dat_3D_converted}



def retrieve_downwelling_in_Wm2(cwv = 10, in_wavelength = False):
    '''
    If in_wavelength is False, dict returned has keys 'wavenumbers' and 'downwelling flux'. Wavenumber in [cm-1], downwelling flux in [W.m-2/cm-1].
    If in_wavelength is True, dict returned has keys 'wavelengths' and 'downwelling flux'. Wavelengths in [um], downwelling flux in [W.m-2/um]

    :param cwv: column water vapor, in mm (either 10, 24 or 54 mm)
    :param in_wavelength: boolean, determines whether downwelling flux is retrieved per wavenumber [cm-1] or per wavelength [um]
    :return: dict with two 1D arrays, as determined by in_wavelength boolean.
    '''
    d_xarray = load_file_as_xarray(cwv)
    x_wavnums = d_xarray['wavenumber']

    if in_wavelength:
        downwelling_flux_cm2_um = convert_from(d_xarray.sel(column='downwelling_flux'),
                                               units_in='per wavenumber [/cm-1]', units_out='per wavelength [/um]',
                                               corresponding_xs = x_wavnums) # [W.cm-2/cm-1] --> [W.cm-2/um]
        downwelling_flux = 1e4 * downwelling_flux_cm2_um  # [W.cm-2/um] --> [W.m-2/um]
        wavns = convert_from(x_wavnums, units_in='wavenumber [cm-1]', units_out='wavelength [um]')
        return {'wavelengths': wavns, 'downwelling flux': downwelling_flux}  # {[um], [W.m-2/um]}

    else:  # in wavenumber
        downwelling_flux = 1e4 * d_xarray.sel(column='downwelling_flux')  # [W.cm-2/cm-1] --> [W.m-2/cm-1]
        wavns = d_xarray['wavenumber']  # [cm-1]

        return {'wavenumbers':wavns, 'downwelling flux':downwelling_flux}  # {[m-1], [W.m-2/m-1]}


def retrieve_downwelling_in_particleflux(cwv=10):
    '''

    :param cwv: column water vapour to retrieve, in mm (either 10, 24 or 54 mm)
    :return: dict containing two 1D arrays, with keys 'photon energies' [in eV] and 'downwelling photon flux' [in W.m-2/eV]
    '''

    # step 1 - convert from wavenumber to eV
    dw_dict = retrieve_downwelling_in_Wm2(cwv, in_wavelength=False)  # [downwelling flux in W.m-2/cm-1]
    phot_energies_eV = convert_from(dw_dict['wavenumbers'], units_in='wavenumber [cm-1]', units_out='photon energy [eV]')
    dwf_W = convert_from(dw_dict['downwelling flux'], units_in='per wavenumber [/cm-1]', units_out='per photon energy [/eV]')

    # step 2 - convert from W to particles per second
    dwf_part = dwf_W / (q*phot_energies_eV)  # [J.s-1.m-2/eV] / [J.eV-1][eV] --> [s-1.m-2/eV]

    return {'photon energies':phot_energies_eV, 'downwelling photon flux':dwf_part}



def planck_dist(Eph, mu, kT):
    '''
    :param Eph: Particular photon energy, in eV
    :param mu: Fermi level splitting, in eV
    :param kT: in eV
    :return: Photon flux [s-1.m-2] for the given eV, calculated using generalized Planck's law (blackbody with temp T if mu=0).
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

def Ndot_downwellheaviside(Eg, downwell_dict):
    Ephs = downwell_dict['photon energies']
    get_heavisided = np.heaviside(Ephs - Eg, 0.5) * downwell_dict['downwelling photon flux']
    pflux = np.trapz(get_heavisided, Ephs)
    return pflux

def neg_powerdensity_downwellheaviside(mu, Eg, Ephs_p, kT_converter, downwell_dict):
    # used to optimize over mu, with scipy
    Ndot_out = Ndot_planckheaviside(Ephs_p, Eg, mu, kT_converter)
    Ndot_in = Ndot_downwellheaviside(Eg, downwell_dict)
    J = q*(Ndot_out-Ndot_in)
    return mu*J

def neg_powerdensity_plancks(mu, Eg, Ephs, kT_convert, kT_env):
    # used to optimize over mu, with scipy
    Ndot_out = Ndot_planckheaviside(Ephs, Eg, mu, kT_convert)
    Ndot_in = Ndot_planckheaviside(Ephs, Eg, 0, kT_env)
    J = q * (Ndot_out - Ndot_in)
    return mu * J


