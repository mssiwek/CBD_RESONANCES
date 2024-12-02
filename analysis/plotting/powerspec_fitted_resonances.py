import numpy as np
import matplotlib.pyplot as plt
import misc
import PowSpec as PS
import seaborn as sns
import pickle as pkl
import utils

figwidth = 10
figheight = 5
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2
Pb = (2.*np.pi)

def plot(all_fps, figpath, titles, param, \
        accs, f_accs, fg_cavs, fname, xlim = [0.1, 20], ylim=[1.e-17, 1.e-7],\
        **kwargs):

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
    
    for ax,fp,title in zip(axes, all_fps, titles):
        # ax.set_title(title, fontsize=fs)
        ax.grid(color='k', alpha=0.1, linestyle='-', linewidth=0.5*linewidth)
        for i,acc in enumerate(accs):
            for j,f_acc in enumerate(f_accs):
                for l,fg_cav in enumerate(fg_cavs):
                    if fg_cav:
                        label = r'$\rm{Cavity \ incl.},$' + '\n'
                        if acc:
                            label += r'$\dot{M}_{\rm b} \geq 0 $'
                        else:
                            label += r'$\dot{M}_{\rm b} = 0 $'
                        if f_acc: 
                            label += r'$, \ f_{\rm acc} \ \rm{incl}$'
                        else:
                            label += r'$, \ f_{\rm acc} \ \rm{excl}$'
                    else:
                        label = r'$\rm{Cavity \ excl.},$' + '\n'
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
                    
                    powspec = PS.PowSpec(t/Pb, val)
                    powspec.ls(xlim=xlim)
                    powspec.ls_peaks(npeaks=0)
                    powspec.nearest_res(lim_m = 10, lim_l = 10, width=0.05)
                    powspec.get_colors()

                    ax.semilogy(np.log10(powspec.x_ls), powspec.y_ls, \
                            color=colors[i][j][l], \
                            label=label, linewidth=linewidth, alpha=0.4)

                    for k,n in enumerate(powspec.match_res['n']):
                        label_res = 'n=%d, %s, l=%d, m=%d,\n log10(f) = %.2f' \
                                                        %(n, powspec.match_res['res_type'][k], \
                                                        powspec.match_res['m'][k],\
                                                        powspec.match_res['l'][k],\
                                                        np.log10(powspec.match_res['res_omega'][k]))
                        
                        ax.vlines(np.log10(powspec.match_res['res_omega'][k]), \
                                    color=powspec.match_res['colors'][k], \
                                    ymin=ylim[0], \
                                    ymax=ylim[1], \
                                    label = label_res, \
                                    linewidth=2*linewidth, alpha=0.4)
        
        ax.set_xlabel(r'$\rm{\log_{10}[frequency}/\omega_{\rm b}]$', fontsize=fs)
        ax.set_ylabel(r'$\rm{Amplitude}$', fontsize=fs)

        ax.set_xlim(np.log10(xlim))
        ax.set_ylim(ylim)
        
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 0.7, 0.5, 0.2), \
                  borderaxespad=0.3, ncol=1, \
                  fancybox=True, shadow=True, fontsize=0.75*fs)
    
    misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth, axes=axes)
    # fig.subplots_adjust(bottom=0.3, wspace=0.33)
    plt.tight_layout()
    plt.savefig(figpath + fname + '.png')
    plt.close()


if __name__ == "__main__":
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
    tmax = 5000*Pb

    #plot ebdot or abdot
    param = 'ebdot'

    # mass accretion included or not
    accs = [False]
    # linear momentum accretion included or not
    f_accs = [False]
    # Include or exclude cavity?
    fg_cavs = [True, False]

    xlim = [0.4, 10]
    ylim = [1.e-17, 1.e-7]

    fname = 'test_powspec'
    if False in fg_cavs:
        fname += '_compare_r>a'
    fname += '_eb=%.2f_qb=%.2f' %(param_dict['e'][0], param_dict['q'][0])


    plot(all_fps, figpath, titles, param, \
        accs, f_accs, fg_cavs, fname, tmin=tmin, tmax=tmax, \
        xlim=xlim, ylim=ylim)