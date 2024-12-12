import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import utils
import plot_utils as put

figwidth = 10
figheight = 7
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2
palette = sns.color_palette('colorblind')
Pb = (2.*np.pi)



""" --------------------------------------------------- """
""" HELPER FUNCTION TO SAVE FIGURES FOR TIMESERIES DATA """
""" --- THIS IS TO HELP AVOID REPEATED CODE SNIPPETS -- """
""" --------------------------------------------------- """
def sfig(fig, figpath, fname):
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    put.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
    plt.tight_layout()
    plt.savefig(figpath + fname + '.png')
    plt.close()
    return()
""" --------------------------------------------------- """
""" HELPER FUNCTION TO SAVE FIGURES FOR TIMESERIES DATA """
""" --- THIS IS TO HELP AVOID REPEATED CODE SNIPPETS -- """
""" --------------------------------------------------- """



""" ---------------------------------------------------- """
""" HELPER FUNCTION TO PLOT TIMESERIES DATA, AND SAVE IT """
""" -- IF NO AXIS IS SUPPLIED OR NEEDS TO BE RETURNED -- """
""" ---------------------------------------------------- """
def plot_ts(fp, label, color, \
            param, acc, f_acc, fg_cav, \
            raw = True, mean=True, smoothed=True, \
            normalized = True, \
            **kwargs):
    
    """ LOAD AXIS """
    if 'ax' in kwargs:
        ax = kwargs['ax']
    else:
        fig,ax = plt.subplots(nrows=1, ncols=1, \
                              sharex=True, sharey=True,\
                              figsize = (figwidth, figheight))
        fname = kwargs['fname']
        figpath = kwargs['figpath']
    


    """ ----------------- """
    """ LOAD ALL THE DATA """
    """ ----------------- """
    binevol = utils.load_binevol(fp + '/BinEvol/', acc, f_acc, fg_cav)

    if param == 'ebdot':
        val = binevol.ebdot()
    if param == 'abdot':
        val = binevol.abdot()

    t = binevol.SimInit.t
    if 'tmin' in kwargs and 'tmax' in kwargs:
        inds_lower = t >= kwargs['tmin']
        inds_upper = t <= kwargs['tmax']
        inds = inds_lower & inds_upper
        t = t[inds]
        val = val[inds]
        mbdot_norm = binevol.mbdot_norm()[inds]
    else:
        t = binevol.SimInit.t
        mbdot_norm = binevol.mbdot_norm()
    
    if normalized:
        val = val/np.mean(mbdot_norm)
        # this is the rate of change of val per unit time,
        # normalized by the mean accretion rate in this time interval
    """ ----------------- """
    """ LOAD ALL THE DATA """
    """ ----------------- """



    if raw:
        """ ------------------------- """
        """ PLOT RAW TIME-SERIES DATA """
        """ ------------------------- """
        ax.plot(t/Pb, \
                val, \
                '-', \
                color=color, \
                label=label, \
                linewidth=linewidth, alpha=0.2)
        """ ------------------------- """
        """ PLOT RAW TIME-SERIES DATA """
        """ ------------------------- """



    if mean:
        """ ----------------------------- """
        """ PLOT MEAN OF TIME-SERIES DATA """
        """ ----------------------------- """
        ax.plot(t/Pb, \
                np.ones(np.shape(val))*np.mean(val),\
                '--',\
                color=color, \
                label = 'mean=%.2e' %(np.mean(val)), \
                linewidth=linewidth, alpha=0.8)
        """ ----------------------------- """
        """ PLOT MEAN OF TIME-SERIES DATA """
        """ ----------------------------- """



    if smoothed:
        """ ------------------------------ """
        """ PLOT SMOOTHED TIME-SERIES DATA """
        """ ------------------------------ """
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
                color=color, \
                linewidth=0.8*linewidth)
        """ ------------------------------ """
        """ PLOT SMOOTHED TIME-SERIES DATA """
        """ ------------------------------ """



    """ ------------------------------ """
    """ FORMAT AXIS, LABELS AND LEGEND """
    """ ------------------------------ """
    ax.grid(color='k', alpha=0.1, linestyle='-', linewidth=0.5*linewidth)
    if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        ax.set_ylim(kwargs['ylim'])
    
    ax.set_yscale('symlog', linthresh=1.e-1)
    
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.2), \
                mode="expand", borderaxespad=0, ncol=2, \
                fancybox=True, shadow=True, fontsize=fs)
    ax.set_xlabel(r'$t/P_{\rm b}$', fontsize=fs)
    if param == 'ebdot':
        ax.set_ylabel(r'$\dot{e}_{\rm b}$', fontsize=fs)
    if param == 'abdot':
        ax.set_ylabel(r'$\dot{a}_{\rm b}$', fontsize=fs)
    """ ------------------------------ """
    """ FORMAT AXIS, LABELS AND LEGEND """
    """ ------------------------------ """

    if 'ax' in kwargs:
        return(ax)
    else:
        sfig(fig,figpath,fname, **kwargs)
        return()
""" ---------------------------------------------------- """
""" HELPER FUNCTION TO PLOT TIMESERIES DATA, AND SAVE IT """
""" -- IF NO AXIS IS SUPPLIED OR NEEDS TO BE RETURNED -- """
""" ---------------------------------------------------- """





""" --------------------------------------------------------- """
""" MAIN PLOTTING ROUTINE FUNCTION THAT IDENTIFIES THE PARAMS """
""" ------- THAT ARE TO BE PLOTTED AGAINST EACH OTHER ------- """
""" --------------------------------------------------------- """
def plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = True, \
        raw = True, mean=True, smoothed=True, \
        figpath = '../figures/', fname='ebdot', **kwargs):
    """Plots the timeseries of param.

    Args:
        all_fps: Filepaths of simulation data that will be plotted.
        labels: list of strings to identify plotted data in legend.
        param: str, either 'ebdot' or 'abdot'
        accs, f_accs, fg_cavs: These are *always* lists which define different 
                               combinations of mass and angular momentum accretion
                               in the simulation,
                               and choose whether to incl/excl the cavity in the 
                               gravitational torque calculation.
        

    Returns:
        nothing.
    """

    """
    
    """

    if all_in_one_ax:
        nrows = 1
        nested = False
    else:
        nrows = len(all_fps)
        if isinstance(all_fps[0], (list,np.ndarray)):
            """ NOW WE HAVE NESTED LISTS OF FPS, 
                AND NEED TO LOOP OVER EACH OF THEM """
            nested = True
        else:
            nested = False

    fig,axes = plt.subplots(nrows=nrows, ncols=1, \
                            sharex=True, sharey=True,\
                            figsize = (figwidth, nrows*figheight))
    if not isinstance(axes, (list,np.ndarray)):
        axes=np.array([axes])

    
    # """ --------------------------------------------- """
    # """ LOAD THE COLORS FOR EACH LINE WE ARE PLOTTING """
    # """ --------------------------------------------- """
    if 'colors' in kwargs:
        colors = kwargs['colors']
    else:
        if all_in_one_ax:
            colors = np.empty([len(all_fps), 3])
            for i,c in zip(np.arange(len(all_fps)),palette):
                colors[i] = c
        elif nested:
            colors = palette
        else:
            """ COMPARING accs, f_accs and fg_cavs """
            #need the 3rd dim = 3 because colors are given in tuples
            colors = np.empty([len(accs), len(f_accs), len(fg_cavs), 3]) 
            k = 0
            for i in np.arange(len(accs)):
                for j in np.arange(len(f_accs)):
                    for l in np.arange(len(fg_cavs)):
                        colors[i][j][l] = palette[k]
                        k+=1
    # """ --------------------------------------------- """
    # """ LOAD THE COLORS FOR EACH LINE WE ARE PLOTTING """
    # """ --------------------------------------------- """



    # """ --------------------------------------------------------------------------------------------- """
    # """ PLOT MULTIPLE RUNS (DIFFERENT EB, QB, PHI, ... ) IN ONE AXIS (WITH ALL OTHER PARAMS THE SAME) """
    # """ --------------------------------------------------------------------------------------------- """
    if all_in_one_ax:
        ax = axes[0]
        acc, f_acc, fg_cav = accs[0], f_accs[0], fg_cavs[0]
        for color,fp,label in zip(colors,all_fps, labels):  
            if 'title' in kwargs:
                ax.text(0.01,0.99, kwargs['title'],
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes, fontsize=fs)                    
            ax = plot_ts(fp, label, color, \
                    param, acc, f_acc, fg_cav, \
                    raw=raw, mean=mean, smoothed=smoothed, \
                    ax=ax, **kwargs)
    # """ --------------------------------------------------------------------------------------------- """
    # """ PLOT MULTIPLE RUNS (DIFFERENT EB, QB, PHI, ... ) IN ONE AXIS (WITH ALL OTHER PARAMS THE SAME) """
    # """ --------------------------------------------------------------------------------------------- """



    # """ ---------------------------------------------------------------------- """
    # """ EACH AXIS PLOTS SEVERAL FPS IN ONE, BUT WITH MULTIPLE AXIS LOOPING OVER
    #     NESTED ARRAYS IN all_fps TO CONTRAST PARAMETERS ACROSS THE SIMULATIONS """
    # """ ---------------------------------------------------------------------- """
    elif nested:
        acc, f_acc, fg_cav = accs[0], f_accs[0], fg_cavs[0]
        for ax,all_fp,title in zip(axes,all_fps,kwargs['titles']):
            for fp,label,color in zip(all_fp,labels,colors):
                ax.text(0.01,0.99, title,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes, fontsize=fs)
                ax = plot_ts(fp, label, color, \
                        param, acc, f_acc, fg_cav, \
                        raw = raw, mean=mean, smoothed=smoothed, \
                        ax = ax, **kwargs)
    # """ ---------------------------------------------------------------------- """
    # """ EACH AXIS PLOTS SEVERAL FPS IN ONE, BUT WITH MULTIPLE AXIS LOOPING OVER
    #     NESTED ARRAYS IN all_fps TO CONTRAST PARAMETERS ACROSS THE SIMULATIONS """
    # """ ---------------------------------------------------------------------- """
    
    
    # """ ----------------------------------------------------------- """
    # """ PLOT (MULTIPLE) RUN(S) (EB, QB, PHI, ... ) IN SEPARATE AXES,
    #     CHOOSING THE COMBINATION OF acc, f_acc and f_cav IN EACH AX """
    # """ ----------------------------------------------------------- """
    else:
        for ax,fp,title in zip(axes,all_fps,kwargs['titles']):
            k = 0
            ax.text(0.01,0.99, title,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax.transAxes, fontsize=fs)
            for i,acc in enumerate(accs):
                for j,f_acc in enumerate(f_accs):
                    for l,fg_cav in enumerate(fg_cavs):
                        label = labels[k]
                        ax = plot_ts(fp, label, colors[i][j][l], \
                                    param, acc, f_acc, fg_cav, \
                                    raw = raw, mean=mean, smoothed=smoothed, \
                                    ax = ax, **kwargs)
                        k+=1

    sfig(fig, figpath, fname)
    return()
# """ --------------------------------------------------------- """
# """ MAIN PLOTTING ROUTINE FUNCTION THAT IDENTIFIES THE PARAMS """
# """ ------- THAT ARE TO BE PLOTTED AGAINST EACH OTHER ------- """
# """ --------------------------------------------------------- """
    

if __name__ == "__main__":
    #cbd_ebqb_study_rs_videos
    fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study/res=5.000e-10/alpha=0.100/dim=2/"
    figpath = "../figures/"

    """ CASE1: PLOTTING MULTIPLE SIMULATIONS IN ONE AXIS """
    CASE1 = True
    """ CASE2: PLOTTING MULTIPLE SIMULATIONS IN ONE AXIS, WITH MULTIPLE SUBPLOTS """
    CASE2 = False
    """ CASE3: PLOTTING ONE SIM PER AXIS, WITH MULTIPLE COMBINATIONS OF 
               f_acc,acc,fg_cav, 
               WITH POTENTIALLY MULTIPLE SUBPLOTS (EACH FOR A DIFFERENT SIMULATION) """
    CASE3 = False


    if CASE1:
        fname = 'CASE1'
        test_param = 'phi'
        param_dict = \
            {\
            'e': [0.3],
            'q': [1.0],
            'aspect_ratio': [0.10],
            'phi': [0.0,180],
            'AccretionRadius': [0.03],
            'SinkRate': [0.50],\
            'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
            }
        all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

        tmin = 2000*Pb
        tmax = 2010*Pb
        ylim = [-500,500]#[-5.e-6, 5.e-6]
        param = 'ebdot'
        # mass accretion included or not
        accs = [False]
        # linear momentum accretion included or not
        f_accs = [False]
        # Include or exclude cavity?
        fg_cavs = [True]

        labels = []
        for tp in param_dict[test_param]:
            if test_param == 'phi':
                labels.append(r'$\Phi = $' + '%d' %(tp))
            if test_param == 'e':
                labels.append(r'$e_{\rm b} = $' + '%.2f' %(tp))
            if test_param == 'q':
                labels.append(r'$q_{\rm b} = $' + '%.2f' %(tp))
        
        title = r'$e_{\rm b} = %.2f, q_{\rm b} = %.2f$' %(param_dict['e'][0], param_dict['q'][0])
        
        plot(all_fps, labels, param, \
            accs, f_accs, fg_cavs, \
            all_in_one_ax = True, \
            raw = True, mean=True, smoothed=True, \
            figpath = figpath, fname=fname, \
            tmin=tmin, tmax=tmax, ylim=ylim, title=title)
    
    if CASE2:
        fname='CASE2'
        all_fps = []
        titles = []
        ebs = [0.1, 0.8]
        test_param = 'q'
        for eb in ebs:
            param_dict = \
                {\
                'e': [eb],
                'q': [0.1,1.0],
                'aspect_ratio': [0.10],
                'phi': [0.0],
                'AccretionRadius': [0.03],
                'SinkRate': [0.50],\
                'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
                }
            all_fp = utils.get_fpaths_params(fp_root, param_dict)
            all_fps.append(all_fp)
            titles.append(r'$e_{\rm b} = %.2f$' %eb)

        tmin = 8000*Pb
        tmax = 8010*Pb
        ylim = [-500,500] #[-5.e-6, 5.e-6]
        param = 'ebdot'
        # mass accretion included or not
        accs = [False]
        # linear momentum accretion included or not
        f_accs = [False]
        # Include or exclude cavity?
        fg_cavs = [True]

        labels = []
        for tp in param_dict[test_param]:
            labels.append(r'$q_{\rm b} = $' + '%.2f' %(tp))

        plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = False, \
        raw = True, mean=True, smoothed=True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim, \
        titles=titles)

    
    if CASE3: 
        fname = 'CASE3'
        test_param = 'e'
        ebs = [0.1,0.8]
        param_dict = \
            {\
            'e': ebs,
            'q': [1.0],
            'aspect_ratio': [0.10],
            'phi': [0.0],
            'AccretionRadius': [0.03],
            'SinkRate': [0.50],\
            'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
            }
        all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

        titles = []
        for eb in ebs:
            titles.append(r'$e_{\rm b} = %.2f,\ q_{\rm b} = %.2f$' %(eb, param_dict['q'][0]))

        tmin = 8000*Pb
        tmax = 8005*Pb
        ylim = [-500,500] #[-5.e-6, 5.e-6]
        param = 'ebdot'
        # mass accretion included or not
        accs = [False]
        # linear momentum accretion included or not
        f_accs = [False]
        # Include or exclude cavity?
        fg_cavs = [True, False]

        labels = []
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
                    labels.append(label)
        
        plot(all_fps, labels, param, \
            accs, f_accs, fg_cavs, \
            all_in_one_ax = False, \
            raw = True, mean=True, smoothed=True, \
            figpath = figpath, fname=fname, \
            tmin=tmin, tmax=tmax, ylim=ylim,\
            titles=titles)