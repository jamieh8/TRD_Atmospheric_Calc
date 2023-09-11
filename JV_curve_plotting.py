import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *

# Drawing I-V, power curves
cmap = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]

Egs_planck = np.arange(0.001, 0.20, 0.02)

cutoff_angle = 90

Egs_dw = np.arange(0.15, 0.20, 0.002)
# Egs = [0.02]

to_plot = []

cwv, Tc = 10, 296.724
emitter = planck_law_body(Tc, Ephs)

# atm_dataset = atmospheric_dataset(cwv)
# combined = TRD_in_atmosphere(emitter, atm_dataset)
# to_plot += [{'label':'cwv10', 'colormap':cmap_pink, 'atmosphere data':atm_dataset, 'TRD':emitter, 'Egs':Egs_dw, 'TRD in atm': combined}]

atm_100K = planck_law_body(100, Ephs)
combined_atm100K = TRD_in_atmosphere(emitter, atm_100K)
to_plot += [{'label':'100K env', 'colormap':cmap_orange, 'atmosphere data':atm_100K, 'TRD':emitter, 'Egs':Egs_planck, 'TRD in atm': combined_atm100K}]


for scheme in to_plot:
    fig, axs = plt.subplots(1,2, layout='tight')
    Egs = scheme['Egs']
    for iE, Eg in enumerate(Egs):
        # fermi level splitting
        mus = np.linspace(-Eg, 0, 100)

        cmap_s = scheme['colormap']
        i_colour = cmap_s(0.2 + 0.8 * iE / len(Egs))

        TRD_obj = scheme['TRD']
        atm_obj = scheme['atmosphere data']
        TRD_in_atm_obj = scheme['TRD in atm']

        # N_out, emitted from converter
        # N_in, emitted from environment, abs by converter
        N_out_mu_planck = np.vectorize(TRD_obj.retrieve_Ndot_heaviside)(Eg, mus, cutoff_angle)
        N_in_mu = atm_obj.retrieve_Ndot_heaviside(Eg, mu=0, cutoff_angle=cutoff_angle)

        # calculate and plot with N abs from env from (1/orange) blackbody at T=30K, (2/pink) downwelling data
        J_mu = q*(N_out_mu_planck - N_in_mu)  # [J] [m-2.s-1] --> A.m-2. J = N_out - N_in

        pds = np.vectorize(TRD_in_atm_obj.power_density)(mus, Eg, cutoff_angle)

        axs[0].plot(mus, J_mu, label = f'E$_g$ = {Eg:.3f} eV', color=i_colour)  # current
        axs[1].plot(mus, pds, color=i_colour)  # power density

    axs[0].set_xlabel("$\mu$ [eV] / V [V]")
    axs[0].set_ylabel('Current, J [A.m$^{-2}$]')
    axs[0].legend()

    axs[1].set_xlabel("$\mu$ [eV] / V [V]")
    axs[1].set_ylabel("Power Density [W.m$^{-2}$]")

    for ax in axs:
        ax.plot([-1, 1], [0, 0], '-k', lw=0.8)
        ax.set_xlim([-Egs[-1], 0])



plt.show()