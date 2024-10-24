import numpy as np
import matplotlib.pyplot as plt
import misc
import pickle as pkl
import BinEvol as BE
import time
import seaborn as sns
import os
import pickle as pkl
import matplotlib 
import scipy
import astropy
from astropy.timeseries import LombScargle

figwidth = 10
figheight = 7
fs = 25
ticksize = 8
tickwidth = 1.5
linewidth = 2
Pb = (2.*np.pi)


def plot(all_fps, figpath, titles, param, \
        accs, f_accs, \
        fname='ebdot', **kwargs):

    start_time_init = time.time()
    fig,axes = plt.subplots(nrows=len(all_fps), ncols=1, \
                                sharex=True, sharey=True,\
                                figsize = (figwidth, len(all_fps)*figheight))
    if not isinstance(axes, (list,np.ndarray)):
        axes=[axes]
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    #need the 3rd dim = 3 because colors are given in tuples
    colors = np.empty([len(accs), len(f_accs), 3]) 
    k = 0
    for i in np.arange(len(accs)):
        for j in np.arange(len(f_accs)):
            colors[i][j] = sns.color_palette('colorblind')[k]
            k+=1

    for fp,ax,title in zip(all_fps, axes, titles):
        # ax.set_title(title, fontsize=fs)
        ax.grid(color='k', alpha=0.1, linestyle='-', linewidth=0.5*linewidth)
        for i,acc in enumerate(accs):
            for j,f_acc in enumerate(f_accs):
                label = ''
                if acc:
                    label += r'$\dot{M}_{\rm b} \geq 0 $'
                else:
                    label += r'$\dot{M}_{\rm b} = 0 $'
                if f_acc: 
                    label += r'$, \ f_{\rm acc} \ \rm{incl}$'
                else:
                    label += r'$, \ f_{\rm acc} \ \rm{excl}$'

                """ CHECK IF OBJECT ALREADY EXISTS """
                binevol_name = 'BinEvoltxt'
                if acc:
                    binevol_name += '_acc'
                else:
                    binevol_name += '_no_acc'
                if f_acc:
                    binevol_name += '_f_acc'
                else:
                    binevol_name += '_no_f_acc'

                filename = fp + '/BinEvol/' + binevol_name + '.pkl'
                if os.path.exists(fp + '/BinEvol/'):
                    if os.path.isfile(filename):
                        print("files exist!")
                        start_time = time.time()
                        filehandler = open(filename, 'rb')
                        print("filehandler = ", filehandler)
                        binevol = pkl.load(filehandler)
                        print("loading %s.pkl took %s seconds ---" %(binevol_name,(time.time() - start_time)))
                    else:
                        print("dir exists, but files do not")
                        siminit = BE.SimInit(fp, acc=acc, f_acc=f_acc)
                        binevol = BE.BinEvol(siminit, from_txt=True, from_snap=False)
                        filehandler = open(filename, 'wb') 
                        pkl.dump(binevol, filehandler)
                else:
                    os.makedirs(fp + '/BinEvol/')
                    start_time = time.time()
                    siminit = BE.SimInit(fp, acc=acc, f_acc=f_acc)
                    print("loading SimInit took %s seconds ---" % (time.time() - start_time))
                    start_time = time.time()
                    binevol = BE.BinEvol(siminit, from_txt=True, from_snap=False)
                    print("loading BinEvol took %s seconds ---" % (time.time() - start_time))
                    filehandler = open(filename, 'wb') 
                    pkl.dump(binevol, filehandler)


                start_time = time.time()
                if param == 'ebdot':
                    val = binevol.ebdot()
                if param == 'abdot':
                    val = binevol.abdot()
                print("calculating param took %s seconds ---" % (time.time() - start_time))

                start_time = time.time()
                if 'tmin' in kwargs and 'tmax' in kwargs:
                    inds_lower = binevol.SimInit.t >= kwargs['tmin']
                    inds_upper = binevol.SimInit.t <= kwargs['tmax']
                    inds = inds_lower & inds_upper
                    t = binevol.SimInit.t[inds]
                    val = val[inds]
                else:
                    t = binevol.SimInit.t

                """
                FOURIER ANALYSIS
                """
                frequency, power = LombScargle(t/Pb, val, normalization='psd').autopower()

                ax.loglog(np.log10(frequency), power, \
                        color=colors[i][j], \
                        label=label, linewidth=linewidth, alpha=0.8)

                # from scipy.signal import savgol_filter
                # # Savitzky-Golay filter
                # # https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                # if 'filter_window' in kwargs:
                #     filter_window = kwargs['filter_window']
                # else:
                #     filter_window = 21
                # val_smoothed = savgol_filter(val, filter_window, 3)
                # ax.plot(t/Pb, \
                #         val_smoothed, \
                #         '-',\
                #         color=colors[i][j], \
                #         linewidth=0.8*linewidth)
        
        ax.set_xlabel(r'$\rm{\log_{10}[frequency}/\omega_{\rm b}]$', fontsize=fs)
        ax.set_ylabel(r'$\rm{Amplitude}$', fontsize=fs)

        start_time = time.time()
        if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs:
            ax.set_ylim(kwargs['ylim'])
        
        # ax.set_yscale('symlog', linthresh=1.e-8)
        
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.2), \
                    mode="expand", borderaxespad=0, ncol=2, \
                  fancybox=True, shadow=True, fontsize=0.7*fs)
        print("Setting limits and plotting legends took %s seconds ---" % (time.time() - start_time))

    misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
    plt.tight_layout()
    start_time = time.time()
    plt.savefig(figpath + fname + '.png')
    print("Saving figure took %s seconds ---" % (time.time() - start_time))
    plt.close()

    print("Entire plotting routine took %s seconds ---" % (time.time() - start_time_init))



            
