from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def fill_in_downwelling(atm, Tskin, low=True):
    Fph_sofar = atm.retrieve_spectral_array(yvals = 's-1.m-2', xvals = 'eV', col_name = 'downwelling_flux')

    Ephs_sofar = atm.photon_energies
    Ephs_after = np.arange(0.31+6.2*1e-5, 1, 6.2*1e-5)

    planck_filler = planck_law_body(T=Tskin)
    Fph_after = planck_filler.spectral_photon_flux(Eph=Ephs_after, mu=0, cutoff_angle=90)

    # from Helen: "moderately" transmissive window from 0.32 - 0.36 eV
    if low:
        correction_fac = 0.95 * np.heaviside(Ephs_after-0.37, 0.5) + 0.05 # is 1 if Eph > 0.37 eV, 0.8 if < 0.37 eV
        Fph_after *= correction_fac

    return {'Ephs':np.append(Ephs_sofar,Ephs_after), 'Fphs':np.append(Fph_sofar, Fph_after)}


Egs = np.linspace(0.05, 0.31, 50)

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)



fig, axs = plt.subplots(1,3, layout='tight')
plot_dict = []

cmap = plt.get_cmap('Set1')
plot_dict = [{'label':'High est (BB fill)', 'fill_str':'high', 'colour':cmap(0)},
             {'label':'Low est (with window)', 'fill_str':'low', 'colour':cmap(1)},
             {'label':'Modelled', 'fill_str':'none', 'colour':cmap(2)}]

atm_dat_sample = atmospheric_dataset_new(cwv='low', location='telfer', Tskin=301.56)
emitter = planck_law_body(T=301.56, Ephs=atm_dat_sample.photon_energies)

for pdct in plot_dict:
    atm_dat = atmospheric_dataset_new(cwv='low', location='telfer', Tskin=301.56, spectral_fill_type=pdct['fill_str'])
    dict = atm_dat.fill_in_downwelling()

    axs[0].plot(dict['Ephs'], dict['Fphs'], label = pdct['label'], c=pdct['colour'])

    Ndots= []
    for Eg in Egs:
        # Fphs_heavisded = dict['Fphs'] * np.heaviside(dict['Ephs']-Eg, 0.5)
        # integral_over_Eph = np.trapz(y=Fphs_heavisded, x=dict['Ephs'])
        integral_over_Eph = atm_dat.retrieve_Ndot_heaviside(Eg, cutoff_angle=None)
        Ndots += [integral_over_Eph]

    axs[1].plot(Egs, Ndots, c=pdct['colour'])

    for out_int, ls in zip(['quad', 'trapz'], ['-', '--']):
        PDs = []
        for Eg in Egs:
            trd_in_env = TRD_in_atmosphere(emitter, atm_dat, out_int_method=out_int)
            opt_pd = trd_in_env.optimize_mu(Eg, cutoff_angle=None)
            PDs += [(-1) * opt_pd['max power']]
        axs[2].plot(Egs, PDs, ls, c=pdct['colour'])


axs[0].set_ylabel('Spectral Photon Flux Density [s$^{-1}$.m$^{-2}$/eV]')
axs[0].set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
axs[0].legend()
axs[0].set_yscale('log')

axs[1].set_ylabel('$\dot{\mathrm{N}}_\mathrm{in}$ [s$^{-1}$.m$^{-2}$]')
axs[1].set_xlabel('Bangap, E$_g$ [eV]')
axs[1].set_yscale('log')

custom_legend = [Line2D([0],[0], color = 'k', linestyle='solid', label='quad '),
                       Line2D([0],[0], color = 'k', linestyle='dashed', label='trapz')]

axs[2].set_ylabel('Max Power Denstiy [W.m$^{-2}$]')
axs[2].set_xlabel('Bangap, E$_g$ [eV]')
axs[2].set_yscale('log')
axs[2].legend(handles=custom_legend, title='emission int. method')

plt.show()