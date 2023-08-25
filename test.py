from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

T = 300
kT = kb * T / q  # kT in eV
print(kT)

dw_dict = retrieve_downwelling_flux(10, in_wavelength=True)  # [downwelling flux in W.m-2/m-1]

# convert from wavelength to eV
phot_energies_J = h * c / (1e-6 * dw_dict['wavelengths'])  # E = hf = hc/lambda,  [m.s-1][J.s]/[m] --> [J]
phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
dwf_W = dw_dict['downwelling flux'] * 1e6 * h * c / (
            phot_energies_J * phot_energies_eV)  # [W.m-2.m-1]*[J.s][m.s-1]*[J-1][eV-1] --> [W.m-2/eV]

# step 1 - convert from wavenumber to eV
# phot_energies_J = h*c*(dw_dict['wavenumbers'])  # E = hf = hc v, [m-1]*[m.s-1][J.s] --> [J]
# phot_energies_eV = phot_energies_J / q  # [J] --> [eV]
# dwf_W = dw_dict['downwelling flux'] * h*c/q  # [W.m-2/m-1]*[][] --> [W.m-2/eV]

# step 2 - convert from W to particles per second
dwf_part = dwf_W / phot_energies_J  # [J.s-1.m-2/eV] --> [s-1.m-2/eV]


ephs = np.arange(0,0.35,0.01)
pd= 2*np.pi/(c**2 * (h/q)**3) * np.vectorize(planck_dist)(ephs, 0, kT) # units [eV^2] --> *q^2 --> [J^2]


# plt.plot(dw_dict['wavelengths'], dw_dict['downwelling flux'])
# plt.xlabel('Wavelength [um]')
# plt.ylabel('Downwelling flux, [W.m-2/um]')


plt.plot(phot_energies_eV, dwf_W)
plt.xlabel('photon E [eV]')
plt.ylabel('W.m-2/eV')


# plt.plot(ephs, pd)
# plt.plot(phot_energies_eV, dwf_part)
# plt.xlabel('photon energies [eV]')
# plt.ylabel('s-1.m-2/eV')

plt.show()