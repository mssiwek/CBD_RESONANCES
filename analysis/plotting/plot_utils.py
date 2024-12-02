import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import misc


figwidth = 7
figheight = 6
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2

def plot_sinks(ax, sn=None, size=0.1, BoxSize=300):
    if sn is not None:
        for sink_no in np.arange(0,len(sn['PartType5']['ParticleIDs'])):
            marker_coords = sn['PartType5']['Coordinates'][sink_no][0:2] - BoxSize/2.
            mtot = np.sum(sn['PartType5']['Masses'])
            msize = size * sn['PartType5']['Masses'][sink_no]/mtot
            circle = plt.Circle((marker_coords[0], marker_coords[1]), msize, \
                                linewidth=1, color='k', fill=True)
            ax.add_patch(circle)
    else:
        print("WARNING: SINK PLOTTING WITHOUT SNAPSHOT NOT YET IMPLEMENTED.")
    
    return(ax)

def plot_mesh(X, Y, Z, **kwargs):
    if 'ax' in kwargs:
        ax = kwargs['ax']
    else:
        fig = plt.figure(figsize = (figwidth, figheight))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        ax = fig.add_subplot(1,1,1)
    
    v_range = 1.e-8
    if 'vlims' in kwargs:
        vlims = kwargs['vlims']
    else:
        vlims = [-1.e-8, 1.e-8]
    linthresh=1.e-10
    if 'linthresh' in kwargs:
        linthresh = kwargs['linthresh']
    
    cbar_norm_log = colors.SymLogNorm(linthresh=linthresh, linscale=0.03, \
                                vmin=vlims[0], vmax=vlims[1])
    
    im = ax.pcolormesh(X, Y, Z, norm=cbar_norm_log, \
                        cmap='RdBu_r')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    try:
        cbar.set_label(kwargs['cbar_label'], rotation=90, fontsize=fs)
    except:
        print("Warning: Did not supply cbar_label in kwargs.")

    ''' Add a marker for the sink particles '''
    plot_sinks(ax)

    ax.set_xlabel(r'$x/a_{\rm b}$', fontsize=fs)
    ax.set_ylabel(r'$y/a_{\rm b}$', fontsize=fs)

    if 'title' in kwargs:
        ax.set_title(kwargs['title'], fontsize=fs)

    plt.tight_layout()

    if 'ax' in kwargs:
        return(ax)
    else:
        misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
        plt.tight_layout()
        if 'figpath' in kwargs and 'fname' in kwargs:
            plt.savefig(kwargs['figpath'] + kwargs['fname'] + '.png')
        else:
            print("No figure path or figure name given! Saving in current path as snap_interp.png")
            plt.savefig('snap_interp.png')
            plt.close()

def lighten_color(color, a):
    """
    a: values over 1 will darken, values under 1 will lighten color
    from: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - a * (1 - c[1]), c[2])
