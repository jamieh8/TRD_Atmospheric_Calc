import numpy as np
import os
import xarray as xr
from solcore.constants import kb, c, h, q

def planck_dist(Eph, mu, kT):
    return Eph**2/(np.exp((Eph-mu)/kT)-1)

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


def retrieve_downwelling_flux(cwv = 10, in_wavelength = False):
    '''

    :param cwv: column water vapor, in mm (either 10, 24 or 54 mm)
    :return: dict with keys 'wavenumber' and 'downwelling flux', values are 1D arrays. Wavenumber in [m-1], downwelling flux in [W.m-2/m-1]
    '''
    d_xarray = load_file_as_xarray(cwv, convert_to_wavelength=in_wavelength)

    if in_wavelength:
        downwelling_flux = 1e-4 * d_xarray.sel(column='downwelling_flux')  # [W.cm-2/um] --> [W.m-2/um]
        wavns = d_xarray['wavelength']  # [um]
        return {'wavelengths': wavns, 'downwelling flux': downwelling_flux}

    else:  # in wavenumber
        downwelling_flux = 1e2 * d_xarray.sel(column='downwelling_flux')  # [W.cm-2/cm-1] --> [W.m-2/m-1]
        wavns = 1e2 * d_xarray['wavenumber']  # [cm-1] --> [m-1]

        return {'wavenumbers':wavns, 'downwelling flux':downwelling_flux}


def retrieve_downwelling_particleflux(cwv=10):
    dw_dict = retrieve_downwelling_flux(cwv, in_wavelength=True)  # [downwelling flux in W.m-2/m-1]

    # convert from wavelength to eV
    phot_energies_J = h*c/(1e-6*dw_dict['wavenumber'])  # E = hf = hc/lambda,  [m.s-1][J.s]/[m] --> [J]
    phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
    dwf_W = dw_dict['downwelling flux'] * 1e6 * h*c/(phot_energies_J*phot_energies_eV)  # [W.m-2.m-1]*[J.s][m.s-1]*[J-1][eV-1] --> [W.m-2/eV]

    # step 1 - convert from wavenumber to eV
    # phot_energies_J = h*c*(dw_dict['wavenumber'])  # E = hf = hc v, [m-1]*[m.s-1][J.s] --> [J]
    # phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
    # dwf_W = dw_dict['downwelling flux'] * h*c/q  # [W.m-2/m-1]*[][] --> [W.m-2/eV]

    # step 2 - convert from W to particles per second
    dwf_part = dwf_W / phot_energies_J  # [J.s-1.m-2/eV] --> [s-1.m-2/eV]

    return [phot_energies_eV, dwf_part]







