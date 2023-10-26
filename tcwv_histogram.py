import numpy as np
import matplotlib.pyplot as plt

copied_data = np.loadtxt('histogram-plot-data.csv', skiprows=1, delimiter=',')

edges = np.arange(0,copied_data[-1,0]+1,1)
print(edges)

fig, axs = plt.subplots(1,2, layout='tight')
axs[0].stairs(copied_data[:,1], edges=edges, fill=True, fc='silver')

vert_line_dicts = [
    {'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'ls': 'dotted', 'symbol': 's'},

                   {'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'ls':'solid', 'symbol':'o'},
                   {'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','ls':'solid', 'symbol':'o'},
                   {'tcwv':70.51, 'Tskin':299.86, 'color':'teal','ls':'solid', 'symbol':'o'},

                   {'tcwv':2.87, 'Tskin':287.31, 'color':'lightblue', 'ls':'dashed', 'symbol':'^'},
                   {'tcwv':19.97, 'Tskin':301.828, 'color':'royalblue', 'ls':'dashed', 'symbol':'^'},
                   {'tcwv':37.91, 'Tskin':299.096, 'color':'darkblue', 'ls':'dashed', 'symbol':'^'}]

for line_dict in vert_line_dicts:
    tcwv = line_dict['tcwv']
    Tskin = line_dict['Tskin']
    insert_index = np.searchsorted(edges, line_dict['tcwv'])
    try:
        yval = copied_data[:,1][insert_index-1]
    except:
        yval = 0
    # axs[0].plot(2*[line_dict['tcwv']], [0,4], ls=line_dict['ls'], c=line_dict['color'])
    axs[0].annotate('', xy=(tcwv, yval), xytext=(tcwv,yval+0.5), arrowprops = dict(arrowstyle='->', color='black'))
    axs[0].plot([tcwv], [yval+0.55], line_dict['symbol'], c=line_dict['color'])

    axs[1].plot([tcwv], [Tskin], line_dict['symbol'], c=line_dict['color'], markersize=7)

ax_scatter = axs[1]
ax_scatter.set_xlabel('TCWV [mm]')
ax_scatter.set_ylabel('Skin Temp [K]')
ax_scatter.minorticks_on()
ax_scatter.grid()

ax_hist = axs[0]
ax_hist.set_xlabel('TCWV [mm]')
ax_hist.set_ylabel('Occurence [%]')
ax_hist.set_ylim([0,4])

plt.show()