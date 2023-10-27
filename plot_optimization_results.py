from TRD_Atmospheric_Functions import *
import matplotlib.pyplot as plt

new_file_dicts = [
    {'loc':'telfer', 'cwvstring':'low', 'tcwv':6.63, 'Tskin':301.56, 'color':'darkorange', 'symbol':'o', 'Eg':0.094},
    {'loc':'telfer', 'cwvstring':'mid', 'tcwv':34.45, 'Tskin':306.43,'color':'darkviolet','symbol':'o', 'Eg':0.094},
    {'loc':'telfer', 'cwvstring':'high', 'tcwv':70.51, 'Tskin':299.86, 'color':'teal','symbol':'o', 'Eg':0.1}
    ]

sweep_str = 'cutoff_angle'
alg_powell = pg.scipy_optimize(method='Powell', tol=1e-5)
alg_de = pg.de(gen=50, ftol=1e-5)

to_plot = []
for AD_dict in new_file_dicts:
    Eg = AD_dict['Eg']
    cwv = AD_dict['tcwv']
    cwv_str = AD_dict['cwvstring']
    loc_str = AD_dict['loc']

    #get file
    filename = f'{loc_str} {cwv_str}_sweep_{sweep_str}.csv'
    data = np.loadtxt(fname = filename, delimiter=',')

    # relative change
    PD_90 = data[:,1][-1]
    y_rel = data[:,1]/PD_90

    # get optimum point
    atm_dat = atmospheric_dataset_new(cwv_str, loc_str)
    TRD = planck_law_body(AD_dict['Tskin'], Ephs = atm_dat.photon_energies)
    combined_trd_env = TRD_in_atmosphere(TRD, atm_dat)
    arg_f = {'Eg':Eg, 'eta_ext':1, 'consider_nonrad':False}
    opt_xs, opt_pd = get_best_pd(combined_trd_env, args_to_opt=['cutoff_angle', 'mu'], args_to_fix=arg_f, alg=alg_powell)
    x_opt = opt_xs['cutoff_angle']
    print(x_opt)
    y_opt = opt_pd[0]/PD_90
    print(opt_pd)
    print(PD_90)

    to_plot += [{'label':f'{loc_str} {cwv_str}', 'color':AD_dict['color'], 'x': data[:,0], 'y': y_rel, 'x opt':x_opt, 'y opt':y_opt}]



for plt_dat in to_plot:
    plt.plot(plt_dat['x'], plt_dat['y'], c=plt_dat['color'])
    plt.plot(plt_dat['x opt'], plt_dat['y opt'], 'o', c=plt_dat['color'])

plt.plot([0,90], [1,1], '--k', lw=1)
plt.xlim([10,90])
# plt.ylim([0,1.1])
plt.minorticks_on()
plt.grid()
plt.ylabel('Relative Power Density, P($\\theta_c$) / P(90$^\circ$)')
plt.xlabel('Cutoff Angle, $\\theta_c$')

plt.show()