import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *


ds_lsts = [get_dataset_list()[0:3], get_dataset_list()[3:6], get_dataset_list()[6:9]]


fig, axs = plt.subplots(2,3, layout='tight')
for ci, ds_lst in enumerate(ds_lsts):
    for ds in ds_lst:
        atm_dat = atmospheric_dataset_new(cwv=ds['cwvstring'], location=ds['loc'], Tskin=ds['Tskin'], date='23dec')
        ys = atm_dat.retrieve_spectral_array(yvals='s-1.m-2', xvals='eV', col_name='downwelling_flux')
        axs[0][ci].plot(atm_dat.photon_energies, ys, c=ds['color'], label=ds['loc']+' '+ds['cwvstring'], alpha=0.8)

        ys = atm_dat.retrieve_spectral_array(yvals='W.m-2', xvals='um', col_name='downwelling_flux')
        axs[1][ci].plot(atm_dat.wavelengths, ys, c=ds['color'], alpha=0.8)

for ax in axs[0][:]:
    ax.set_yscale('log')
    ax.set_xlabel('Photon energy, $E_g$ [eV]')
    ax.set_ylabel('Spectral Photon Flux Density, F$_{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
    ax.legend()
    ax.grid()

for ax in axs[1][:]:
    ax.set_xlim([0,30])
    ax.set_xlabel('Wavelength [um]')
    ax.set_ylabel('Spectral Irradiance, F$_e$ [W.m$^{-2}$/um]')
    ax.grid()



plt.show()



