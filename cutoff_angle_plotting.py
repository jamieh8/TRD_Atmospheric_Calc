import numpy as np
import matplotlib.pyplot as plt
from TRD_Atmospheric_Functions import *

new_file_dicts = get_dataset_list()[0:3]

for case in new_file_dicts:
    loc, tcwvstr = case['loc'], case['cwvstring']
    filename = f'{loc} {tcwvstr}_sweep_cutoff_angle.csv'
    loaddat = np.loadtxt(fname = filename, delimiter=',')
    angles, powers = loaddat[:,0], loaddat[:,1]
    relative_change = powers/powers[-1]
    plt.plot(angles,relative_change, c=case['color'], label=f'{loc.capitalize()} {tcwvstr}', lw=2)

plt.plot([0,90], [1,1], '--k', lw=1)
plt.xlim([10,90])
plt.ylim([0,1.1])

plt.xlabel('Cutoff Angle, $\\theta_c$ [°]')
plt.ylabel('Relative Power, P($\\theta_c$) / P(90°)')

plt.legend(loc='upper left')
plt.minorticks_on()
plt.grid()

plt.show()
