import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar


# Tc = 300  # temperature of emitter / converter
mu = 0
kT_eV = kb / q

Te = 30 # temperature of environment
kT_e_eV = kb * Te / q
k_eV = kb/q



Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]



col_dict = {'planck_emit':'skyblue', 'planckheavy_emit':'darkblue', 'planck_env_30K':'orangered', 'planck_env_100K':'orange',
            'boltzmann_emit':'springgreen', 'downwell_env':'lightpink', 'downwellheavy_env':'deeppink'}


# Comparing 3 cwvs
# Optimize over bandgap for full-hemisphere emission and absorption
# plt.figure(1)
# fig3, axs3 = plt.subplots(1,2, layout='tight')
# Egs_dw = np.arange(0.062, 0.20, 0.002)
# cwvs = [10,24,54]
# cwv_cols = {10: 'deeppink', 24:'steelblue', 54:'seagreen'}
# cwv_temps = {10:296.724, 24:304.868, 54:303.512}
#
# for cwv in cwvs:
#     atm_data = atmospheric_dataset(cwv)
#     downwelling_diff = atm_data.retrieve_spectral_array('s-1.m-2','eV','downwelling_flux')
#     # plt.figure(1)
#     # plt.plot(pflux_downwell['photon energies'], pflux_downwell['downwelling photon flux'], label=f'cwv {cwv}', c=cwv_cols[cwv])
#
#     maxPs = []
#     Vmpps = []
#     for Eg in Egs_dw:
#         kT_c_eV = cwv_temps[cwv] * kT_eV
#         opt_mu_dwh = minimize_scalar(neg_powerdensity_downwellheaviside, bounds=[-Eg, 0],
#                                      args=(Eg, Ephs, kT_c_eV, atm_data.photon_energies, downwelling_diff))
#         maxPs += [opt_mu_dwh.fun]
#         Vmpps += [opt_mu_dwh.x]
#
#     axs3[0].plot(Egs_dw, maxPs, label=f'{cwv_temps[cwv]}K out, cwv {cwv} in', c=cwv_cols[cwv])
#     axs3[1].plot(Egs_dw, Vmpps, label=f'{cwv_temps[cwv]}K out, cwv {cwv} in', c=cwv_cols[cwv])


# plt.legend()
# plt.ylabel('Photons density flux [s$^{-1}$.m$^{-2}$/eV')
# plt.xlabel('Photon energy, E$_{ph}$ [eV]')



# Plotting max power density vs Eg - cwv=10 compared to 30K, 100K blackbody environment
# fig3, axs3 = plt.subplots(1,2, layout='tight')
# # neg_powerdensity_downwellheaviside(mu, Eg, Ephs, kT_c_eV, pflux_downwell)
# # find best mu for given Eg
# Egs_dw = np.arange(0.062, 0.20, 0.002)
# Egs_planck = np.arange(0, 0.20, 0.002)
# pflux_downwell = retrieve_downwelling_in_particleflux(10)
#
# maxPs_dwh = []
# maxPs_30K = []
# maxPs_100K = []
# Vmpps = {'dwh':[], '30K':[], '100K':[]}
#
# for Eg in Egs_dw:
#     opt_mu_dwh = minimize_scalar(neg_powerdensity_downwellheaviside, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, pflux_downwell))
#     maxPs_dwh += [opt_mu_dwh.fun]
#     Vmpps['dwh']+= [opt_mu_dwh.x]
#
# for Eg in Egs_planck:
#     opt_mu_30K = minimize_scalar(neg_powerdensity_plancks, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, k_eV*30))
#     maxPs_30K += [opt_mu_30K.fun]
#     Vmpps['30K'] += [opt_mu_30K.x]
#
#     opt_mu_100K = minimize_scalar(neg_powerdensity_plancks, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, k_eV * 100))
#     maxPs_100K += [opt_mu_100K.fun]
#     Vmpps['100K'] += [opt_mu_100K.x]
#
# axs3[0].plot(Egs_dw, maxPs_dwh, color=col_dict['downwellheavy_env'], label='300K out, Atmospheric+Heaviside in')
# axs3[0].plot(Egs_planck, maxPs_30K, color=col_dict['planck_env_30K'], label='300K out, 30K+Heaviside in')
# axs3[0].plot(Egs_planck, maxPs_100K, color='orange', label='300K out, 100K+Heaviside in')
#
# axs3[1].plot(Egs_dw, Vmpps['dwh'], label='Atmospheric DW data', color = col_dict['downwellheavy_env'])
# for k, c_str in zip(['30K', '100K'], ['planck_env_30K', 'planck_env_100K']):
#     axs3[1].plot(Egs_planck, Vmpps[k], color = col_dict[c_str])


# axs3[0].set_ylabel('Power Density [W.m$^{-2}$]')
# axs3[0].set_xlabel('Bandgap, E$_g$ [eV]')
# axs3[0].set_title('Max Power Density vs Bandgap')
# axs3[0].legend()
#
# axs3[1].plot([-1,1],2*[-1*kT_c_eV], '--k')
# axs3[1].text(x=0.150, y=-kT_c_eV+0.0005, s='-kT$_c$')
# axs3[1].set_xlim([0.06,0.2])
# axs3[1].set_xlabel('Bangap, E$_g$ [eV]')
# axs3[1].set_ylabel('V$_{mpp}$ [V]')
# axs3[1].set_title('Optimal Voltage (V$_{mpp}$ vs Bandgap)')


# check results from optimizer, at a particular Eg, cwv
# plt.figure()
# Eg = 0.075
# pf_dw_cwv10 = retrieve_downwelling_in_particleflux(10)
# mus = np.arange(-Eg, 10e-5, 0.001)
# power_vs_mu = []
# for mu in mus:
#     power_vs_mu += [neg_powerdensity_downwellheaviside(mu, Eg, Ephs, kT_c_eV, pf_dw_cwv10)]
# plt.plot(mus, power_vs_mu)
#
# opt_mu_dwh = minimize_scalar(neg_powerdensity_downwellheaviside, bounds=[-Eg, 0], args=(Eg, Ephs, kT_c_eV, pf_dw_cwv10))
# plt.plot(opt_mu_dwh.x, opt_mu_dwh.fun, 'o')


# Test: optimize over cutoff angle
cwv, T = 10, 296.724

atm_data = atmospheric_dataset(cwv)
downwelling_diff = atm_data.retrieve_spectral_array('s-1.m-2', 'eV', 'downwelling_flux')
emitter_planck = planck_law_body(T=T)

combined_trd_env = TRD_in_atmosphere(emitter_planck, atm_data, Ephs)
alg = pg.scipy_optimize(method='SLSQP', tol=1e-5)
# alg = pg.de(gen=50, ftol=1e-3)

max_pds = []
for angle in np.arange(10,95,10):
    opt_res = get_best_pd(combined_trd_env, cutoff_angle=angle, alg = alg)
    max_pds += [opt_res[1]]
    print(f'Eg {opt_res[0][0]}, mu {opt_res[0][1]}')

plt.plot(np.arange(10,95,10), max_pds)


plt.show()