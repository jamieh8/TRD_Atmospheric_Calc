from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator


# ----- Import EQE and escape probability arrays -----

# load measured, temperature dependent EQE from csv
csv_name = 'PVI-5_9427_T_EQE' #'PVI-6_8967_T_EQE'#'PVIA-5_14032_T_EQE' #"PVI-4_7868_T_EQE.csv"  # 'PVIA-10_14065_T_EQE.csv'
fn = fr"C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Temp dependent EQE\{csv_name}.csv"

# first row is temperature, first column is wavelength. values are EQE in %.
data = np.genfromtxt(fn, dtype=float, delimiter=',')
wavelengths = data[1:,0] # in [nm]
Ts = data[0,1:]  # [K]
EQEs = data[1:,1:]/100

# set up interpolator so that EQE can be retrieved at arbirarity temp
EQE_interpolator = RegularGridInterpolator((wavelengths, Ts), EQEs, bounds_error=False, fill_value=None)

# --- (uncomment to plot the 2D EQE) ---
# fig, hEQE_ax = plt.subplots(1,1,layout='tight')
# hmap = hEQE_ax.pcolormesh(Ts, wavelengths*1e-3, data[1:,1:], cmap='magma', shading='nearest')
# hEQE_ax.set_xlabel('Diode temp [K]')
# hEQE_ax.set_ylabel('Wavelength [um]')
# cbar = plt.colorbar(hmap)
# cbar.ax.set_ylabel('EQE [%]')

# load ray tracing results as 2 arrays (modelled emission angles beta, escape_prob(beta))
fn = r"C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Andreas Dropbox - Calculations\hyperhemi_T_withinangle_u_45.txt"
raytrace_data = np.loadtxt(fn, dtype=float)
betas_mod = raytrace_data[:,0]
escape_prob = raytrace_data[:,1]


# ----- Jsc calc & plotting -----

# set up diode and environment temps to sweep
K_spacing = 1
diode_temps = np.arange(min(Ts), max(Ts), K_spacing)
env_temps = np.arange(100, 300, K_spacing)

# convert units, [nm]->[um]->[eV]
# not required, but the integration happens over Eph so doing this once, up front, speeds things up slightly
Ephs = convert_from(wavelengths*1e-3, units_in='wavelength [um]', units_out='photon energy [eV]')

fn_results = f'Jsc_{csv_name}_{K_spacing}Kspacing.csv'
filepath_results = r'C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Temp dependent EQE\\' + fn_results
plot_from_file_only = True


if plot_from_file_only:
    # don't calculate, just load in results
    load_dat = np.loadtxt(filepath_results, skiprows=1, delimiter=',')
    currents = load_dat[:,1:]
else:
    # calculate results
    currents = [], []
    for diode_T in diode_temps:
        # retrieve interpolated diode EQE at this temperature
        diodeEQE_atT = EQE_interpolator((wavelengths, diode_T))
        np.nan_to_num(diodeEQE_atT, copy=False)

        # set up diode_EQE object using interpolated EQE and escape probability (dynamic resistance could be added here for power calcs)
        diode = diode_EQE(spec_type='measured', spec_xs={'spec x':Ephs, 'spec units':'photon energy [eV]'}, spec_EQE=diodeEQE_atT,
                                    angle_type='escape probability - isotropic', escape_prob=escape_prob, betas=betas_mod)

        emitter_b = mathemetica_style_body(T=diode_T)  # define emitter temperature (emission profile)

        row=[]
        for env_T in env_temps:
            env_b = mathemetica_style_body(T=env_T)  # define environment temperature
            die = diode_in_environment(emitter=emitter_b, environment=env_b, diode_EQE = diode)  # combine
            Jsc = die.current_density(V=0)  # retrieve current density at V=0 (s.c.)
            row += [Jsc]

        currents += [row]

    # save results
    ar_currents =  np.array(currents)
    data_to_save = np.hstack((np.transpose(np.array(diode_temps, ndmin=2)), ar_currents))  # add diode T as first col
    data_to_save = np.vstack((np.insert(env_temps, 0, 0), data_to_save))  # add 0+env temp as first row
    np.savetxt(fn_results, data_to_save, delimiter=',')


# plotting Jsc heatmap
fig, h_ax = plt.subplots(1,1,layout='tight')
hmap = h_ax.pcolormesh(env_temps, diode_temps, currents, cmap='magma', shading='nearest')
cbar = plt.colorbar(hmap)
h_ax.set_ylabel('Diode Temperature [K]')
h_ax.set_xlabel('Environment Temperature [K]')
cbar.ax.set_ylabel('J$_{sc}$ [A.m$^{-2}$]')

h_ax.set_title(csv_name)


plt.show()
