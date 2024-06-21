from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


Egs = np.linspace(0.05, 0.3, 50)

alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

dataset = get_dataset_list()[0]
atm_dat_sample = atmospheric_dataset_new(cwv=dataset['cwvstring'], location=dataset['loc'], Tskin=dataset['Tskin'], date='23dec')
emitter = planck_law_body(T=dataset['Tskin'], Ephs=atm_dat_sample.photon_energies)
ref_atm_high = atmospheric_dataset_new(cwv=dataset['cwvstring'], location=dataset['loc'], Tskin=dataset['Tskin'], spectral_fill_type='high', date='23dec')
trd_in_env_ref = TRD_in_atmosphere(emitter, ref_atm_high, out_int_method='quad')

plot_dict = []

cmap = plt.get_cmap('Set1')
plot_dict = [{'label':'Highest estimate (BB fill)', 'fill_str':'high', 'colour':cmap(0)},
             # {'label':'Mid estimate (with window)', 'fill_str':'low', 'colour':cmap(1)},
             {'label':'Lowest estimate (modelled only)', 'fill_str':'none', 'colour':cmap(2)}]

# cmap = plt.get_cmap('plasma')
# max_Ephs = np.append([0.31],np.arange(0.4, 1, 0.1))
# for i, max_Eph in enumerate(max_Ephs):
#     plot_dict += [{'label':f'to {max_Eph:.2f} eV', 'fill_str':'high', 'max_Eph_fill':max_Eph, 'colour':cmap(i/len(max_Ephs))}]

PDs = []
for Eg in Egs:
    opt_pd = trd_in_env_ref.optimize_mu(Eg, cutoff_angle=None)
    PDs += [opt_pd['max power']]
PDs_ref = (-1)*np.array(PDs)

fig, axs = plt.subplots(1,3, layout='tight')


# for max_Eph in max_Ephs:
for pdct in plot_dict:
    atm_dat = atmospheric_dataset_new(cwv=dataset['cwvstring'], location=dataset['loc'], Tskin=dataset['Tskin'],
                                      spectral_fill_type=pdct['fill_str'], date='23dec')

    dct = atm_dat.fill_in_downwelling()

    axs[0].plot(dct['Ephs'], dct['Fphs'], label = pdct['label'], c=pdct['colour'])

    for out_int, ls in zip(['quad'], ['-']):
        PDs, Vs = [], []
        for Eg in Egs:
            trd_in_env = TRD_in_atmosphere(emitter, atm_dat, out_int_method=out_int)
            opt_pd = trd_in_env.optimize_mu(Eg, cutoff_angle=None)
            PDs += [(-1) * opt_pd['max power']]
            Vs += [opt_pd['Vmpp']]
        axs[1].plot(Egs, PDs, ls, c=pdct['colour'])

    # diff_in_PD = PDs/PDs_ref #100*(PDs-PDs_ref)/PDs_ref
    diff_in_PD = PDs - PDs_ref
    axs[2].plot(Egs, diff_in_PD, c=pdct['colour'])
    # axs[2].text(x=Egs[-1], y=diff_in_PD[-1], s=' ' + pdct['label'], color=pdct['colour'], ha='right', va='bottom')


axs[0].set_ylabel('Spectral Photon Flux Density [s$^{-1}$.m$^{-2}$/eV]')
axs[0].set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
axs[0].legend()
axs[0].set_yscale('log')
add_wn_ticks(axs[0])

axs[1].set_ylabel('Max Power Density [W.m$^{-2}$]')
axs[1].set_xlabel('Bangap, E$_g$ [eV]')
axs[1].set_yscale('log')
add_wn_ticks(axs[1])

axs[2].set_ylabel('Power Density Difference [W.m$^{-2}$]')
axs[2].set_xlabel('Bangap, E$_g$ [eV]')
axs[2].set_yscale('log')
axs[2].set_xlim([0.05,0.36])


plt.show()
