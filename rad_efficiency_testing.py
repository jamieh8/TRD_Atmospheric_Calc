import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *


# Drawing I-V, power curves
cmap_plasma = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']
cmap_purple = matplotlib.colormaps['Purples']

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

Egs_planck = np.arange(0.001, 0.20, 0.02)

cutoff_angle = 90

Egs_dw = np.arange(0.15, 0.20, 0.002)
# Egs = [0.02]

to_plot = []

cwv, Tc = 10, 296.724
emitter = planck_law_body(Tc, Ephs)

atm_dataset = atmospheric_dataset(cwv)
combined_cwv10 = TRD_in_atmosphere(emitter, atm_dataset)
# to_plot += [{'label':'cwv10, $\eta_{rad}$=1', 'color':'darkgreen', 'Eg':0.1, 'TRD in atm': combined_cwv10,
#              'nonrad':True, 'rad efficiency':1}]
# to_plot += [{'label':'cwv10, $\eta_{rad}$=0.5', 'color':'limegreen', 'Eg':0.1, 'TRD in atm': combined_cwv10,
#              'nonrad':True, 'rad efficiency':0.5}]

atm_3K = planck_law_body(3, Ephs)
combined_atm3K = TRD_in_atmosphere(emitter, atm_3K)
rad_effs = np.arange(0.1, 1.1, 0.15)
print(rad_effs)
for ir, rad_eff in enumerate(rad_effs):
    col = cmap_purple(0.2 + 0.8 * ir / len(rad_effs))
    to_plot += [{'label':'3K env, $\eta_{rad}$'+ f'={rad_eff:.2f}', 'color':col, 'Eg':0.1, 'TRD in atm': combined_atm3K,
                 'nonrad':True, 'rad efficiency':rad_eff}]


fig, axs = plt.subplots(1,2, layout='tight')
for scheme in to_plot:
    Eg = scheme['Eg']
    nr_bool = scheme['nonrad']
    rad_eff = scheme['rad efficiency']

    # fermi level splitting
    mus = np.linspace(-Eg, 0, 100)

    i_colour = scheme['color']
    TRD_in_atm_obj = scheme['TRD in atm']

    Js_mu = np.vectorize(TRD_in_atm_obj.current_density)(mus, Eg, cutoff_angle, eta_ext=rad_eff, consider_nonrad = nr_bool)
    pds = Js_mu*mus
    # pds = np.vectorize(TRD_in_atm_obj.power_density)(mus, Eg, cutoff_angle, eta_ext=rad_eff, consider_nonrad = nr_bool)

    axs[0].plot(mus, Js_mu, label = scheme['label'], color=i_colour)  # current
    axs[1].plot(mus, pds, color=i_colour)  # power density

    axs[0].set_xlabel('$\mu$ [eV] / V [V]')
    axs[0].set_ylabel('Current, J [A.m$^{-2}$]')
    axs[0].legend()

    axs[1].set_xlabel('$\mu$ [eV] / V [V]')
    axs[1].set_ylabel('Power Density [W.m$^{-2}$]')

    for ax in axs:
        ax.plot([-1, 1], [0, 0], '-k', lw=0.8)
        ax.set_xlim([-Eg, 0])

axs[0].set_ylim([0,1500])
axs[1].set_ylim([-20,0])

plt.show()