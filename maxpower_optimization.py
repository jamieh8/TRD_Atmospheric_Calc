import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *
from scipy.optimize import minimize_scalar


Tc = 300  # temperature of emitter / converter
mu = 0
kT_c_eV = kb * Tc / q

Te = 30 # temperature of environment
kT_e_eV = kb * Te / q
k_eV = kb/q



Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]



col_dict = {'planck_emit':'skyblue', 'planckheavy_emit':'darkblue', 'planck_env_30K':'orangered', 'planck_env_100K':'orange',
            'boltzmann_emit':'springgreen', 'downwell_env':'lightpink', 'downwellheavy_env':'deeppink'}



# Comparing 3 cwvs
plt.figure(1)
fig3, axs3 = plt.subplots(1,2, layout='tight')
Egs_dw = np.arange(0.062, 0.20, 0.002)
cwvs = [10,24,54]
# maxPs = {'cwv10':[], 'cwv24':[], 'cwv54':[]}
# Vmpps = {'cwv10':[], 'cwv24':[], 'cwv54':[]}
cwv_cols = {10: 'deeppink', 24:'steelblue', 54:'seagreen'}

for cwv in cwvs:
    pflux_downwell = retrieve_downwelling_in_particleflux(cwv)
    plt.figure(1)
    plt.plot(pflux_downwell['photon energies'], pflux_downwell['downwelling photon flux'], label=f'cwv {cwv}', c=cwv_cols[cwv])

    maxPs = []
    Vmpps = []
    for Eg in Egs_dw:
        opt_mu_dwh = minimize_scalar(neg_powerdensity_downwellheaviside, bounds=[-Eg, 0],
                                     args=(Eg, Ephs, kT_c_eV, pflux_downwell))
        maxPs += [opt_mu_dwh.fun]
        Vmpps += [opt_mu_dwh.x]

    axs3[0].plot(Egs_dw, maxPs, label=f'300K out, cwv {cwv} in', c=cwv_cols[cwv])
    axs3[1].plot(Egs_dw, Vmpps, label=f'300K out, cwv {cwv} in', c=cwv_cols[cwv])


plt.legend()
plt.ylabel('Photons density flux [s$^{-1}$.m$^{-2}$/eV')
plt.xlabel('Photon energy, E$_{ph}$ [eV]')



# Plotting max power density vs Eg
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



axs3[0].set_ylabel('Power Density [W.m$^{-2}$]')
axs3[0].set_xlabel('Bandgap, E$_g$ [eV]')
axs3[0].set_title('Max Power Density vs Bandgap')
axs3[0].legend()

axs3[1].plot([-1,1],2*[-1*kT_c_eV], '--k')
axs3[1].text(x=0.150, y=-kT_c_eV+0.0005, s='-kT$_c$')
axs3[1].set_xlim([0,0.2])
axs3[1].set_xlabel('Bangap, E$_g$ [eV]')
axs3[1].set_ylabel('V$_{mpp}$ [V]')
axs3[1].set_title('Optimal Voltage (V$_{mpp}$ vs Bandgap)')

plt.show()