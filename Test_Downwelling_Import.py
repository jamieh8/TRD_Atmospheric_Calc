import numpy as np

from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

def planck_spectral_rad(wav, kT_J):
    wav = wav*1e-6  # [um] --> [m]
    return 2*h*c**2 / (wav**5) * (np.exp(h*c/(wav*kT_J)-1))**(-1)  # [J.s][m.s-1]^2 [m]^-5 --> [J][m2][s-1][m-5] --> [J][m-3][s-1] --> [W.m-3] * per solid angle!!

T = 297
kT_J = kb * T
kT_eV = kT_J / q  # kT in eV [J / (J/eV)] --> [eV]


# Wavelength VS W.m-2/um
dw_dict = retrieve_downwelling_in_Wm2(10, in_wavelength=True)

wvs = dw_dict['wavelengths']
plt.plot(dw_dict['wavelengths'], dw_dict['downwelling flux'], label='downwelling data')
planck_y = 1e-6 * np.vectorize(planck_spectral_rad)(wvs, kT_J) # [W.m-3] --> [W.m-2/um]
plt.plot(wvs, planck_y, label=f'planck, T={T}K')
plt.xlabel('Wavelength [um]')
plt.ylabel('(Spectral) Power Density [W.m-2/um]')
plt.legend()
plt.show()



partflux_dict = retrieve_downwelling_in_particleflux(10)

ephs = np.arange(0,0.35,0.01)  # in [eV]
pd= 2*np.pi/(c**2 * (h/q)**3) * np.vectorize(planck_dist)(ephs, 0, kT_eV)
# units [m.s-1]-2[J.s/J.eV-1]^-3[eV^2] --> [m-2.s2][eV-3.s-3][eV^2] --> [m-2.s-1.eV-1]

plt.plot(partflux_dict['photon energies'], partflux_dict['downwelling photon flux'], label='downwelling data')
plt.plot(ephs, pd, label=f'blackbody, T={T}K')
plt.xlabel('photon energies [eV]')
plt.ylabel('Photon flux [s$^{-1}$.m$^{-2}$/eV]')
plt.legend()
plt.show()

