import numpy as np
import matplotlib.pyplot as plt
import misc
import time
import seaborn as sns
import utils

figwidth = 10
figheight = 7
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2
Pb = (2.*np.pi)


def plot(all_fps, figpath, titles, param, \
        accs, f_accs, fg_cavs, \
        fname='ebdot', **kwargs):

    start_time_init = time.time()
    fig,axes = plt.subplots(nrows=len(all_fps), ncols=1, \
                                sharex=True, sharey=True,\
                                figsize = (figwidth, len(all_fps)*figheight))
    if not isinstance(axes, (list,np.ndarray)):
        axes=[axes]
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    #need the 3rd dim = 3 because colors are given in tuples
    colors = np.empty([len(accs), len(f_accs), len(fg_cavs), 3]) 
    k = 0
    for i in np.arange(len(accs)):
        for j in np.arange(len(f_accs)):
            for l in np.arange(len(fg_cavs)):
                colors[i][j][l] = sns.color_palette('colorblind')[k]
                k+=1

    for fp,ax,title in zip(all_fps, axes, titles):
        # ax.set_title(title, fontsize=fs)
        ax.grid(color='k', alpha=0.1, linestyle='-', linewidth=0.5*linewidth)
        for i,acc in enumerate(accs):
            for j,f_acc in enumerate(f_accs):
                for l,fg_cav in enumerate(fg_cavs):
                    label = ''
                    if fg_cav:
                        if acc:
                            label += r'$\dot{M}_{\rm b} \geq 0 $'
                        else:
                            label += r'$\dot{M}_{\rm b} = 0 $'
                        if f_acc: 
                            label += r'$, \ f_{\rm acc} \ \rm{incl}$'
                        else:
                            label += r'$, \ f_{\rm acc} \ \rm{excl}$'
                    else:
                        label += r'$ f_{\rm g, \, r>a}$'

                    binevol = utils.load_binevol(fp + '/BinEvol/', acc, f_acc, fg_cav)

                    if param == 'ebdot':
                        val = binevol.ebdot()
                    if param == 'abdot':
                        val = binevol.abdot()

                    if 'tmin' in kwargs and 'tmax' in kwargs:
                        inds_lower = binevol.SimInit.t >= kwargs['tmin']
                        inds_upper = binevol.SimInit.t <= kwargs['tmax']
                        inds = inds_lower & inds_upper
                        t = binevol.SimInit.t[inds]
                        val = val[inds]
                    else:
                        t = binevol.SimInit.t
                    ax.plot(t/Pb, \
                            val, \
                            '-', \
                            color=colors[i][j][l], \
                            label=label, linewidth=linewidth, alpha=0.2)
                    ax.plot(t/Pb, \
                            np.ones(np.shape(val))*np.mean(val),\
                            '--',\
                            color=colors[i][j][l], \
                            label = 'mean=%.2e' %(np.mean(val)), linewidth=linewidth, alpha=0.8)

                    from scipy.signal import savgol_filter
                    # Savitzky-Golay filter
                    # https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                    if 'filter_window' in kwargs:
                        filter_window = kwargs['filter_window']
                    else:
                        filter_window = 21
                    val_smoothed = savgol_filter(val, filter_window, 3)
                    ax.plot(t/Pb, \
                            val_smoothed, \
                            '-',\
                            color=colors[i][j][l], \
                            linewidth=0.8*linewidth)
        
        start_time = time.time()
        if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs:
            ax.set_ylim(kwargs['ylim'])
        
        ax.set_yscale('symlog', linthresh=1.e-8)
        
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.2), \
                    mode="expand", borderaxespad=0, ncol=2, \
                  fancybox=True, shadow=True, fontsize=0.7*fs)
        ax.set_xlabel(r'$t/P_{\rm b}$', fontsize=fs)
        print("Setting limits and plotting legends took %s seconds ---" % (time.time() - start_time))

    plt.tight_layout()
    misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
    start_time = time.time()
    plt.savefig(figpath + fname + '.png')
    print("Saving figure took %s seconds ---" % (time.time() - start_time))
    plt.close()

    print("Entire plotting routine took %s seconds ---" % (time.time() - start_time_init))


if __name__ == "__main__":
    #cbd_ebqb_study_rs_videos
    fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study/res=5.000e-10/alpha=0.100/dim=2/"
    figpath = "../figures/"

    param_dict = \
        {\
        'e': [0.2],
        'q': [1.0],
        'aspect_ratio': [0.10],
        'phi': [0.0],
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    
    all_fps = misc.get_fpaths_params(fp_root, param_dict, test_param = 'e')
    titles = []
    for eb in param_dict['e']:
        titles.append(r'$e_{\rm b} = $' + '%.2f' %eb)

    Pb = 2.*np.pi
    tmin = 3000*Pb
    tmax = 3002*Pb

    #note: give xlim in units of Pb
    ylim = [-5.e-6, 5.e-6]#[-5.e-6, 5.e-6]

    #plot ebdot or abdot
    param = 'ebdot'

    # mass accretion included or not
    accs = [False]
    # linear momentum accretion included or not
    f_accs = [False]
    # Include or exclude cavity?
    fg_cavs = [True, False]

    fname = 'timeseries'
    if fg_cavs[0] is False and len(fg_cavs) == 1:
        fname += '_r>a'
    fname += '_%s_qb%.2f' %(param,param_dict['q'][0])
    for eb in param_dict['e']:
        fname += '_%.2feb' %eb
    
    fname = fname + "_accrad=%.3f" %(param_dict['AccretionRadius'][0]) \
            + '_%d<t<%d' %(tmin/Pb, tmax/Pb)

    plot(all_fps, figpath, titles, param, \
             accs, f_accs, fg_cavs, \
             fname=fname, ylim=ylim, \
             tmin=tmin, tmax=tmax)