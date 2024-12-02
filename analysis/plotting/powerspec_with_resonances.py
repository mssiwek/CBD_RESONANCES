import numpy as np
import matplotlib.pyplot as plt
import misc
import pickle as pkl
import BinEvol as BE
import time
import seaborn as sns
import os
import pickle as pkl
from astropy.timeseries import LombScargle
import resonances as res

figwidth = 10
figheight = 7
fs = 25
ticksize = 8
tickwidth = 1.5
linewidth = 2
Pb = (2.*np.pi)

def log10_omega(r):
    return(np.log10(res.omega(r)))

def r_from_omega_log10(omlog, mb=1, G=1):
    om = 10**(omlog)
    return((G*mb/(om**2))**(1/3))


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
                labels_res = ['CR', 'OLR', 'ILR']

                colors_res = ['r', 'b', 'g']
                # for j in enumerate(labels_res):
                #     colors_res.append(sns.color_palette('colorblind')[j])

                r_i = 0
                # for m in range(0,3):
                #     for l in [m-1, m, m+1, m+2, m+3]:

                # ls = [1,2]
                # ms = 3*np.array(ls)
                # for m,l in zip(ms,ls):
                ms = range(0,3)
                mls = [4,5,6,7]
                alphas = np.linspace(0.7,0.1,len(ms)*len(mls))
                alpha_i = 0
                for m in ms:
                    for ml in mls:
                        alpha = alphas[alpha_i]
                        alpha_i += 1
                        l = ml*m
                        if l == 0:
                            continue
                        else:
                            if m > 0:
                                CR = res.omega(res.r_CR(l,m))
                            else:
                                CR = None
                            OLR = res.omega(res.r_OLR(l,m))
                            if m > 1:
                                ILR = res.omega(res.r_ILR(l,m))
                            else:
                                ILR = None
                            
                            RES = [CR, OLR, ILR]

                            for R,label_res,color in zip(RES, labels_res, colors_res):
                                if R is not None:
                                    R_log = np.log10(R)
                                    if (R_log >= kwargs['xlim'][0] and R_log <= kwargs['xlim'][1]):
                                        label_res_i = '%s, m=%d, l=%d' %(label_res, m,l)
                                        # label_res_i += '\n R = %.2f' %(r_from_omega_log10(R_log, mb=1, G=1))
                                        if R_log == 0.00:
                                            continue 
                                        elif abs(R_log) < 0.001:
                                            print("%s, l=%d, m=%d, at log10(R) = %.2f" %(label_res, l, m, R_log))
                                        else:
                                            print("%s, l=%d, m=%d, at log10(R) = %.2f" %(label_res, l, m, R_log))
                                            ax.vlines(R_log, ymin=kwargs['ylim'][0], \
                                                    ymax=kwargs['ylim'][1], \
                                                    color = color, \
                                                    label = label_res_i, \
                                                    linewidth=2*linewidth, \
                                                    alpha=alpha)
                                            r_i += 1

                frequency, power = LombScargle(t/Pb, val, normalization='psd').autopower()

                ax.semilogy(np.log10(frequency), power, \
                        color=colors[i][j], \
                        label=label, linewidth=linewidth, alpha=0.8)

        
        ax_r = ax.secondary_xaxis("top", functions=(r_from_omega_log10, log10_omega))
        from matplotlib.ticker import AutoMinorLocator
        ax_r.xaxis.set_minor_locator(AutoMinorLocator())
        ax_r.set_xlabel(r'$R_{m,l}/a_{\rm b}$', fontsize=fs)
        ax_r.tick_params(labelsize = fs)
        axes.append(ax_r)

        ax.set_xlabel(r'$\rm{\log_{10}[frequency}/\omega_{\rm b}]$', fontsize=fs)
        ax.set_ylabel(r'$\rm{Amplitude}$', fontsize=fs)

        start_time = time.time()
        if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs:
            ax.set_ylim(kwargs['ylim'])
        
        # ax.set_yscale('symlog', linthresh=1.e-8)
        
        # ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.2), \
        #             mode="expand", borderaxespad=0, ncol=2, \
        #           fancybox=True, shadow=True, fontsize=0.7*fs)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.1, 1.0, 0.2), \
                  borderaxespad=0.3, ncol=1, \
                  fancybox=True, shadow=True, fontsize=0.5*fs)
        print("Setting limits and plotting legends took %s seconds ---" % (time.time() - start_time))

    misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth, axes=axes)
    # fig.subplots_adjust(bottom=0.3, wspace=0.33)
    plt.tight_layout()
    start_time = time.time()
    plt.savefig(figpath + fname + '.png')
    print("Saving figure took %s seconds ---" % (time.time() - start_time))
    plt.close()

    print("Entire plotting routine took %s seconds ---" % (time.time() - start_time_init))



if __name__ == "__main__":
    param_dict = \
    {\
    'e': [0.6],
    'q': [1.0],
    'aspect_ratio': [0.10],
    'phi': [0.0],
    'AccretionRadius': [0.03],
    'SinkRate': [0.50],\
    'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
    }

    param = 'ebdot'

    #cbd_ebqb_study_rs_videos
    fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study/res=5.000e-10/alpha=0.100/dim=2/"
    figpath = "../figures/"

    all_fps = misc.get_fpaths_params(fp_root, param_dict, test_param = 'e')
    titles = []
    for eb in param_dict['e']:
        titles.append(r'$e_{\rm b} = $' + '%.2f' %eb)

    tmin = 2000*Pb
    tmax = 5000*Pb

    fname = 'powerspec_%s_qb%.2f' %(param,param_dict['q'][0])
    for eb in param_dict['e']:
        fname += '_%.2feb' %eb 

    fname = fname + "_accrad=%.3f" %(param_dict['AccretionRadius'][0]) \
            + '_%d<t<%d' %(np.ceil(tmin/Pb), np.ceil(tmax/Pb))


    # mass accretion included or not
    accs = [True]
    # linear momentum accretion included or not
    f_accs = [True]

    xlim = [-3,1.5]
    ylim = [1.e-17, 1.e-7]

    plot(all_fps, figpath, titles, param, \
            accs, f_accs, \
            fname=fname,  \
            tmin=tmin, tmax=tmax, \
            xlim=xlim, ylim=ylim)
            
