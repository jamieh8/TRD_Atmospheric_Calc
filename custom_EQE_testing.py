from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt


# filename = os.path.join('simulations_telfer', f'telfer_australia_cwv{cwv}.txt')
filename = r"C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Temp dependent EQE\EQE_5um_diode.csv"
data = np.loadtxt(filename, skiprows=3, delimiter=',')

wavelengths = data[:,0]
EQE = data[:,3]/100

emitter = planck_law_body()
# # wavelengths = np.arange(2, 100)
# Lph_per_um = emitter.angle_spectral_photon_flux(wavelengths, 0, spectral_units='wavelength [um]')
# Fph_per_um = emitter.spectral_photon_flux(wavelengths, 0, spectral_units='wavelength [um]')
#
# # plt.plot(wavelengths, Lph_per_um, label=r'$L_\text{ph Planck}$')
# plt.plot(wavelengths, Lph_per_um*EQE/100, label = r'$L_\text{ph emit}(\lambda, \theta=0)$ (no angular dep)', c='red')
# plt.plot(wavelengths, Fph_per_um*EQE/100, label = r'$F_\text{ph emit}(\lambda)$, Lambertian over $\theta$', c='purple')
# # plt.plot(wavelengths, EQE)
# plt.xlabel('Wavelength, $\lambda$ [um]')
# plt.ylabel('(directional) spectral photon flux [$\mathrm{s^{-1}.m^{-2}.-/um}$]')
# plt.legend()


test_EQE = diode_EQE(spec_type='measured', spec_xs={'spec x':wavelengths, 'spec units':'wavelength [um]'}, spec_EQE=EQE)
tE_xy = test_EQE.get_spec_EQE()
plt.plot(tE_xy['x'], tE_xy['EQE(x)'], label= 'PV-5 measured EQE')

hv_EQE =  diode_EQE(spec_type='heaviside', spec_xs={'spec x':wavelengths, 'spec units':'wavelength [um]'}, Eg = 0.25)
hvE_xy = hv_EQE.get_spec_EQE()
plt.plot(hvE_xy['x'], hvE_xy['EQE(x)'], label ='Heaviside EQE')

print(emitter.retrieve_Ndot(EQE=test_EQE))
print(emitter.retrieve_Ndot(EQE=hv_EQE))
print(emitter.retrieve_Ndot_heaviside(Eg=0.25))


plt.xlabel('Wavelength, $\lambda$ [um]')
plt.ylabel('EQE [%]')
plt.legend()
plt.show()