from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

def fill_in_downwelling(atm, Tskin, low=True):
    Fph_sofar = atm.retrieve_spectral_array(yvals = 's-1.m-2', xvals = 'eV', col_name = 'downwelling_flux')

    Ephs_sofar = atm.photon_energies
    Ephs_after = np.arange(0.31+6.2*1e-5, 1, 6.2*1e-5)

    planck_filler = planck_law_body(T=Tskin)
    Fph_after = planck_filler.spectral_photon_flux(Eph=Ephs_after, mu=0, cutoff_angle=90)

    if low:
        

    return {'Ephs':np.append(Ephs_sofar,Ephs_after), 'Fphs':np.append(Fph_sofar, Fph_after)}



alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

atm_dat = atmospheric_dataset_new(cwv='low', location='telfer')

emitter = planck_law_body(T=300, Ephs = atm_dat.photon_energies)
env = planck_law_body(T=3)


spec_out = fill_in_downwelling(atm_dat, Tskin=301.56, low=False)
plt.plot(spec_out['Ephs'], spec_out['Fphs'])

# Egs = np.linspace(0, 0.3)
# Ndots_trapz  = []
# Ndots_quad = []
# Ndots_env = []
# for Eg in Egs:
#     Ndt = emitter.retrieve_Ndot_heaviside(Eg, mu=0, int_method='trapz')
#     Ndots_trapz += [Ndt]
#
#     Ndq = emitter.retrieve_Ndot_heaviside(Eg, mu=0, int_method='quad')
#     Ndots_quad += [Ndq]
#
#     Nde = atm_dat.retrieve_Ndot_heaviside(Eg, cutoff_angle=None)
#     Ndots_env += [Nde]
#
# plt.plot(Egs, Ndots_trapz, label='emitter, trapz')
# plt.plot(Egs, Ndots_quad, label='emitter, quad')
# plt.plot(Egs, Ndots_env, label='LBLRTM env, trapz')
#
# plt.legend()
#
# plt.ylabel('$\dot{N}$ [s$^{-1}$.m$^{-2}$]')
# plt.xlabel('Bandgap, E$_g$ [eV]')



plt.yscale('log')

plt.show()