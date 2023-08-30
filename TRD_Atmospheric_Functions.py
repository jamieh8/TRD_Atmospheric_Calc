import numpy as np
import os
import xarray as xr
from solcore.constants import kb, c, h, q


def load_file_as_xarray(cwv = 10, convert_to_wavelength = True):
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

    if convert_to_wavelength:
        # converting wavenumber to wavelength: just 1/wavenumber. But also need to adjust flux/radiance since these
        # are per cm-1 (Jacobian transformation)
        ds = 1e-4 * ds * (ds['wavenumber']) ** 2  # [W.cm-2/cm-1]*[cm-1]^2 = [W.cm-2/cm] --> *1e-4 --> [W.cm-2/um]
        ds = ds.rename({'wavenumber': 'wavelength'})
        ds.coords['wavelength'] = 1e4 / ds.coords['wavelength']  # [1/cm-1]=[cm] ---> *1e4 --> [um]

    return ds


def retrieve_downwelling_in_Wm2(cwv = 10, in_wavelength = False):
    '''
    If in_wavelength is False, dict returned has keys 'wavenumbers' and 'downwelling flux'. Wavenumber in [m-1], downwelling flux in [W.m-2/m-1].
    If in_wavelength is True, dict returned has keys 'wavelengths' and 'downwelling flux'. Wavelengths in [um], downwelling flux in [W.m-2/um]

    :param cwv: column water vapor, in mm (either 10, 24 or 54 mm)
    :param in_wavelength: boolean, determines whether downwelling flux is retrieved per wavenumber [m-1] or per wavelength [um]
    :return: dict with two 1D arrays, as determined by in_wavelength boolean.
    '''
    d_xarray = load_file_as_xarray(cwv, convert_to_wavelength=in_wavelength)

    if in_wavelength:
        downwelling_flux = 1e4 * d_xarray.sel(column='downwelling_flux')  # [W.cm-2/um] --> [W.m-2/um]
        wavns = d_xarray['wavelength']  # [um]
        return {'wavelengths': wavns, 'downwelling flux': downwelling_flux}  # {[um], [W.m-2/um]}

    else:  # in wavenumber
        downwelling_flux = 1e2 * d_xarray.sel(column='downwelling_flux')  # [W.cm-2/cm-1] --> [W.m-2/m-1]
        wavns = 1e2 * d_xarray['wavenumber']  # [cm-1] --> [m-1]

        return {'wavenumbers':wavns, 'downwelling flux':downwelling_flux}  # {[m-1], [W.m-2/m-1]}


def retrieve_downwelling_in_particleflux(cwv=10):
    '''

    :param cwv: column water vapour to retrieve, in mm (either 10, 24 or 54 mm)
    :return: dict containing two 1D arrays, with keys 'photon energies' [in eV] and 'downwelling photon flux' [in W.m-2/eV]
    '''

    # step 1 - convert from wavelength to eV
    # dw_dict = retrieve_downwelling_flux(10, in_wavelength=True)  # [downwelling flux in W.m-2/um]
    # phot_energies_J = h * c / (1e-6 * dw_dict['wavelengths'])  # E = hf = hc/lambda,  [m.s-1][J.s]/[m] --> [J]
    # phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
    # dwf_W = dw_dict['downwelling flux'] * 1e6 * h * c / (
    #             phot_energies_J * phot_energies_eV)  # [W.m-2.m-1]*[J.s][m.s-1]*[J-1][eV-1] --> [W.m-2/eV]

    # step 1 - convert from wavenumber to eV
    dw_dict = retrieve_downwelling_in_Wm2(cwv, in_wavelength=False)  # [downwelling flux in W.m-2/m-1]
    phot_energies_J = h * c * (dw_dict['wavenumbers'])  # E = hf = hc v, [m-1]*[m.s-1][J.s] --> [J]
    phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
    dwf_W = dw_dict['downwelling flux'] / (
                h * c / q)  # [W.m-2/m-1]*[J.s][m.s-1][eV/J] --> [W.m-2/m-1]/[eV.m] --> [W.m-2/eV]

    # step 2 - convert from W to particles per second
    dwf_part = dwf_W / phot_energies_J  # [J.s-1.m-2/eV] --> [s-1.m-2/eV]

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


