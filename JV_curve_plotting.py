import matplotlib.pyplot as plt
import matplotlib
from TRD_Atmospheric_Functions import *

# Drawing I-V, power curves

fig2, axs2 = plt.subplots(2,2, layout='tight')
cmap = matplotlib.colormaps['plasma']
cmap_orange = matplotlib.colormaps['Oranges']
cmap_pink = matplotlib.colormaps['RdPu']

Ephs = np.arange(1e-6, 0.31, 0.0001)  # [eV]
Egs_dw = np.arange(0.062, 0.20, 0.002)
Egs_planck = np.arange(0.001, 0.20, 0.002)

Tc = 300  # temperature of emitter / converter
mu = 0
kT_c_eV = kb * Tc / q

Te = 30 # temperature of environment
kT_e_eV = kb * Te / q
k_eV = kb/q

cwv=10
atm_dataset = atmospheric_dataset(cwv)
phot_flux_atm = atm_dataset.retrieve_spectral_array('s-1.m-2', 'eV', 'downwelling_flux')

for ci, Egs in enumerate([Egs_planck, Egs_dw]):
    for iE, Eg in enumerate(Egs[::10]):
        # fermi level splitting
        mus = np.linspace(-Eg, 0, 100)

        # N_out, emitted from converter
        # N_in, emitted from environment, abs by converter
        N_out_mu_planck = np.vectorize(Ndot_planckheaviside, excluded=['Ephs'])(Ephs= Ephs, Eg=Eg, mu=mus, kT=kT_c_eV)
        if ci==0:
            # Planck
            N_in_mu = Ndot_planckheaviside(Ephs, Eg, 0, kT_e_eV)
            i_colour = cmap_orange(0.2 + 0.8 * iE / len(Egs[::10]))
        elif ci==1:
            # Downwelling modelling data
            N_in_mu = Ndot_downwellheaviside(Eg, atm_dataset.photon_energies, phot_flux_atm)
            i_colour = cmap_pink(0.2 + 0.8 * iE / len(Egs[::10]))

        # calculate and plot with N abs from env from (1/orange) blackbody at T=30K, (2/pink) downwelling data
        J_mu = q*(N_out_mu_planck - N_in_mu)  # [J] [m-2.s-1] --> A.m-2. J = N_out - N_in

        axs2[ci][0].plot(mus, J_mu, label = f'E$_g$ = {Eg:.3f} eV', color=i_colour)  # current
        axs2[ci][1].plot(mus, mus*J_mu, color=i_colour)  # power density


for ri in [0,1]:
    for ci in [0,1]:
        axs2[ri][ci].plot([-1,1], [0,0], '-k', lw=0.8)
        axs2[ri][ci].set_xlim([-Egs[-1], 0])


    axs2[ri][0].set_xlabel("Fermi level splitting, $\mu$ [eV] / V [V]")
    axs2[ri][0].set_ylabel("Current, J [A.m$^{-2}$]")
    axs2[ri][0].set_title('Current vs $\mu$')
    axs2[ri][0].legend()

    axs2[ri][1].set_xlabel("Fermi level splitting, $\mu$ [eV] / V [V]")
    axs2[ri][1].set_ylabel("Power Density [W.m$^{-2}$]")
    axs2[ri][1].set_title('Power Density vs $\mu$')


plt.show()