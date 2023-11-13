import numpy as np
import matplotlib.pyplot as plt

new_file_dicts = [
    {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o', 'Eg':0.094},
    {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o', 'Eg':0.094},
    {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o', 'Eg':0.101},

    # {'loc':'california', 'cwvstring':'low', 'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'mid', 'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'symbol': 's'},
    # {'loc':'california', 'cwvstring':'high', 'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'symbol': 's'},
    #
    # {'loc':'tamanrasset', 'cwvstring':'low', 'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'mid', 'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'symbol':'^'},
    # {'loc':'tamanrasset', 'cwvstring':'high', 'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'symbol':'^'}
    ]

for case in new_file_dicts:
    loc, tcwvstr = case['loc'], case['cwvstring']
    filename = f'{loc} {tcwvstr}_sweep_cutoff_angle.csv'
    loaddat = np.loadtxt(fname = filename, delimiter=',')
    angles, powers = loaddat[:,0], loaddat[:,1]
    relative_change = powers/powers[-1]
    plt.plot(angles,relative_change, c=case['color'], label=f'{loc.capitalize()} {tcwvstr}')

plt.plot([0,90], [1,1], '--k', lw=1)
plt.xlim([10,90])
plt.ylim([0,1.1])

plt.xlabel('Cutoff Angle, $\\theta_c \;[\circ]$')
plt.ylabel('Relative Power, P($\\theta_c$) / P(90$^\circ$)')

plt.legend(loc='upper left')
plt.minorticks_on()
plt.grid()

plt.show()
