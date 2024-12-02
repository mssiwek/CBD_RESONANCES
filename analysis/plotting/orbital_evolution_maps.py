import h5py as h5
import gadget
import numpy as np
import matplotlib.pyplot as plt
import sys
import misc
# setting path
sys.path.append('../')
sys.path.append('/n/home00/msiwek/MBHB_disks/arepo_postprocessing/')
import accretion as ac
import accretion as ac
#import json
try:
   import cPickle as pkl
except:
   import pickle as pkl
import matplotlib as mpl
import matplotlib.colors as colors

#plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import accretion as ac

figwidth = 7
figheight = 6
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2


def plot(all_fps, figpath, titles, param, box=20, BoxSize=300, \
    dim='2D', SemiMajorAxis=1, fname='ebdot', contour=False, ratio = False, **kwargs):
    center=[BoxSize/2., BoxSize/2.]

    i_initial, i_final = misc.snap_inds(all_fps, kwargs)


    all_params = []
    for fp in all_fps:
        all_param_fp = open(fp + '/../' + 'all_param.pkl', 'rb')
        all_param = pkl.load(all_param_fp)
        all_param['BoxSize'] = BoxSize
        all_param['acc_txt'] = ac.read_accretion_txt(fp)
        all_params.append(all_param)

    for i in range(i_initial, i_final+1):
        ''' Make figure and axes '''
        if 'ax' not in kwargs:
            n_ax = len(all_fps)
            fig = plt.figure(figsize = (n_ax*figwidth, figheight))
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        sn_time = gadget.Simulation(all_fps[0] + '/snap_%03d.hdf5' %i)
        all_param['snap_time'] = sn_time.Time.value

        vals = {}

        for k,(fp,all_param) in enumerate(zip(all_fps,all_params)):
            ''' First load snapshots '''
            sn = h5.File(fp +  '/snap_%03d.hdf5' %i, 'r')
            sn_sim = gadget.Simulation(fp + '/snap_%03d.hdf5' %i)
            if 'ax' in kwargs:
                ax = kwargs['ax']
            else:
                ax = fig.add_subplot(1, n_ax, k+1)

            ''' Calculate param for each gas cell in domain '''
            vals['%s_%d' %(param,k)], *_ = eval('cean.calc_%s(all_param, sn)' %param)
            ''' Calculate the sum of ebdot and abdot '''
            vals['%s_sum_%d' %(param,k)] = eval('cean.calc_sum_%s(all_param, sn)' %param)

            ax.set_title(titles[k] + r'$, t = %.2f\, P_{\rm b}$' %(sn_time.Time.value/(2*np.pi)) + "\n" + \
                        r'$\dot{%s}_{\rm b} = %.2f \dot{M}_{\rm b}/M_{\rm b}$' %(param[0], vals['%s_sum_%d' %(param,k)]), \
                        fontsize=fs)

            v_range = 1.e-2
            cbar_norm_log = colors.SymLogNorm(linthresh=1.e-5, linscale=0.03, \
                                    vmin=-v_range, vmax=v_range)
            cbar_norm_lin = mpl.colors.Normalize(vmin=-v_range, vmax=v_range)
            
            im = sn_sim.plot_Aslice(vals['%s_%d' %(param,k)], center = center, \
                                    contour=contour, colorbar = False, axes = ax, box = box, \
                                    norm=cbar_norm_log, \
                                    vmin=-v_range, vmax=v_range, \
                                    cmap='RdBu_r')
            
            ax.set_xlabel(r'$x/a_{\rm b}$')
            ax.set_ylabel(r'$y/a_{\rm b}$')

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label(r'$\dot{%s}_{\rm b} \, [\dot{M}_{\rm b}/M_{\rm b}]$' \
                           %(param[0]), rotation=90, fontsize=fs)

            ''' Add a marker for the sink particles '''
            for sink_no in np.arange(0,len(sn['PartType5']['ParticleIDs'])):
                marker_coords = sn['PartType5']['Coordinates'][sink_no][0:2]
                circle = plt.Circle((marker_coords[0], marker_coords[1]), 0.05, linewidth=3, color='r', fill=False)
                ax.add_patch(circle)
                #ax.Circle((marker_coords[0], marker_coords[1]), 0.05, linewidth=3, color='r', fill=False)
                #ax.plot(marker_coords[0], marker_coords[1], 'o', markersize = 0.5, color='r')
            plt.tight_layout()

            if 'ax' in kwargs:
                return(ax)

        if 'ax' not in kwargs:
            misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
            plt.savefig(figpath + fname + '_%03d.png' %i)
            plt.close()

    # misc.make_video(figpath, fname)