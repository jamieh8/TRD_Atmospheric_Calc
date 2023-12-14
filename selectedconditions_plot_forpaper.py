from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def retrieve_dat(month):
    folder = r'C:\Users\z5426944\OneDrive - UNSW\Documents\Thermoradiative Diode\Mean tcwv data from Helen'
    filename = f'{month}_0000UTC_mean.dat'

    dat_array = np.loadtxt(os.path.join(folder, filename))
    dat_array = dat_array.flatten()
    dat_array = dat_array.reshape((721,1440))
    dat_array_spliced = np.append(dat_array[:,720:], dat_array[:,0:720], axis=1)

    return dat_array_spliced

def add_location_annots(ax):
    annots = [{'label':'Tamanrasset', 'coords':(5.5,22.75), 'offset':(-10,30), 'textcoords':(8,55),'symbol':'^', 'textalign':{'ha':'left', 'va':'center'}},
              {'label':'California', 'coords':(-120,34), 'offset':(0, 30), 'textcoords':(-102,64),'symbol':'s', 'textalign':{'ha':'left', 'va':'center'}},
              {'label':'Telfer', 'coords':(122.25, -21.75), 'offset':(-20,-20), 'textcoords':(90,-50),'symbol':'o', 'textalign':{'ha':'right', 'va':'center'}}]
    for annot in annots:
        offset_coords = (annot['coords'][0]+annot['offset'][0], annot['coords'][1]+annot['offset'][1])
        ax.annotate(' ', xy=annot['coords'], xytext=offset_coords, ha='center', va='center', color='black', arrowprops={'arrowstyle':'->', 'color':'red', 'linewidth':2})
        ax.plot(offset_coords[0], offset_coords[1], marker=annot['symbol'], c='white', markersize=15, mec='red', mew=2)
        ax.text(x=annot['textcoords'][0],y=annot['textcoords'][1], s=annot['label'], **annot['textalign'], bbox = {'color':'white', 'pad':0.4})


datasets = get_dataset_list()[0:3]

fig = plt.figure(figsize=(15,7))
gsmap = fig.add_gridspec(2, 3, width_ratios=[3,2,4], hspace=0.1, wspace=0.3)
# gs2 = fig.add_gridspec(1, 3, width_ratios=[5,1,2], wspace=0.05)

ax_mapjan = fig.add_subplot(gsmap[0,0], projection=ccrs.PlateCarree())
ax_mapjune = fig.add_subplot(gsmap[1,0], projection=ccrs.PlateCarree())
map_axs = [ax_mapjan, ax_mapjune]



# TCWV map
jan_dat = retrieve_dat('jan')
june_dat = retrieve_dat('jun')

longs = np.arange(-180, 180, 0.25)
lats = np.arange(90, -90.25, -0.25)

contours = np.arange(0,66+6,6)

for i, dat in enumerate([jan_dat, june_dat]):
    map_ax = map_axs[i]
    map_ax.coastlines()

    hmap = map_ax.contourf(longs, lats, dat,  levels=contours, cmap='rainbow')

    map_ax.set_xticks(np.arange(-180, 180.1, 20))
    map_ax.set_yticks(np.arange(-80, 80.1, 20))
    map_ax.grid(c='white', lw=0.2)

    map_ax.set_xticklabels([])
    map_ax.set_yticklabels([])
    map_ax.xaxis.set_ticks_position('none')
    map_ax.yaxis.set_ticks_position('none')

    map_ax.set_aspect('auto')

add_location_annots(map_axs[0])
map_axs[0].text(x=-180+5,y=-90+3,s='JAN 0000 UTC', ha='left', va='bottom', color='white', bbox = {'color':'black', 'pad':1})
map_axs[1].text(x=-180+5,y=-90+3,s='JUNE 0000 UTC', ha='left', va='bottom', color='white', bbox = {'color':'black', 'pad':1})

fig.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.9)
cb_ax = fig.add_axes([0.05, 0.1, 0.25, 0.02])  # x0, y0, width, height

cbar = fig.colorbar(hmap, cax=cb_ax, orientation='horizontal')
cbar.ax.set_xlabel('Decadal mean TCWV [mm]', fontsize=10)
cbar.ax.set_xticks(contours)
cbar.ax.set_xlim([0, 66])


# add rest of subplots
gs1 = fig.add_gridspec(2, 3, width_ratios=[3,2,4], wspace=0.3, hspace=0.5)
ax_hist = fig.add_subplot(gs1[0,1])
ax_scatter = fig.add_subplot(gs1[1,1])
ax_dwf = fig.add_subplot(gs1[:,2])
axs = [ax_hist, ax_scatter, ax_dwf]

# Histogram and scatter plot
copied_data = np.loadtxt('histogram-plot-data.csv', skiprows=1, delimiter=',')
edges = np.arange(0,copied_data[-1,0]+1,1)
ax_hist.stairs(copied_data[:,1], edges=edges, fill=True, fc='silver')

vert_line_dicts = [
    {'tcwv': 5.32, 'Tskin': 276.298, 'color': 'pink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 17.21, 'Tskin': 295.68, 'color': 'hotpink', 'ls': 'dotted', 'symbol': 's'},
    {'tcwv': 40.32, 'Tskin': 299.231, 'color': 'crimson', 'ls': 'dotted', 'symbol': 's'},

    {'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'ls':'solid', 'symbol':'o'},
    {'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','ls':'solid', 'symbol':'o'},
    {'tcwv':70.51, 'Tskin':299.86, 'color':'mediumseagreen','ls':'solid', 'symbol':'o'},

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
    ax_hist.annotate('', xy=(tcwv, yval), xytext=(tcwv,yval+0.5), arrowprops = dict(arrowstyle='->', color='black'))
    ax_hist.plot([tcwv], [yval+0.55], line_dict['symbol'], c=line_dict['color'])

    ax_scatter.plot([tcwv], [Tskin], line_dict['symbol'], c=line_dict['color'], markersize=7)

ax_scatter.set_xlabel('TCWV [mm]')
ax_scatter.set_ylabel('Skin Temp [K]')
ax_scatter.set_xlim([0,75])
ax_scatter.set_ylim([275, 310])
ax_scatter.minorticks_on()
ax_scatter.grid()

ax_hist.set_xlabel('TCWV [mm]')
ax_hist.set_ylabel('Occurence [%]')
ax_hist.set_ylim([0,4])
ax_hist.set_xlim([0,75])



# Downwelling flux
for ds in datasets:
    cwv_str = ds['cwvstring']
    loc_str = ds['loc']
    atm_data = atmospheric_dataset_new(cwv=cwv_str, location=loc_str, Tskin=ds['Tskin'])
    Ephs = atm_data.photon_energies
    downwelling_ph_fl = atm_data.retrieve_spectral_array(yvals='W.cm-2', xvals='cm-1')
    ax_dwf.plot(Ephs, downwelling_ph_fl, c=ds['color'], label=f'{loc_str.capitalize()} {cwv_str}')

ax_dwf.set_ylabel('Spectral Photon Flux Density, F$_\mathrm{ph}$ [s$^{-1}$.m$^{-2}$/eV]')
# ax.set_ylabel('Spectral Irradiance, F$_e$ [W.cm$^{-2}$/cm$^{-1}$]')
ax_dwf.set_xlabel('Photon Energy, E$_\mathrm{ph}$ [eV]')
ax_dwf.legend(loc='lower left')
# ax.set_xlim([0,0.31])
# ax.set_ylim([0,5*1e23])
ax_dwf.set_yscale('log')

add_wl_ticks(ax_dwf)
# add_wn_ticks(ax)

map_axs[0]. set_title('a)', loc='left', fontsize=15, fontweight='bold', y=1.05, x=-0.1)
for letter, ax in zip(['b','c','d'], axs):
    ax.set_title(letter+')', loc='left', fontsize=15, fontweight='bold', y=1.05, x=-0.1)

plt.show()
