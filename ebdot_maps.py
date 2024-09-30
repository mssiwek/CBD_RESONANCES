import h5py as h5
import glob
import gadget
from gadget import calcGrid
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
import misc
# setting path
sys.path.append('../')
sys.path.append('/n/home00/msiwek/MBHB_disks/arepo_postprocessing/')
import accretion as ac
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
import math
import momentum as mt
import accretion as ac
#import json
try:
   import cPickle as pkl
except:
   import pickle as pkl
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.colors as colors


#plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import viscosity as visc
import pandas as pd
import torques as tq
import accretion_theory as act
import torques_eb_lb as tq_el
import binary_evolution as be
import accretion as ac
import disk_eccentricity as de
import torques as tq
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import inspector_gadget_sliceplots as igs

figwidth = 10
figheight = 6
fs = 30
ticksize = 10
tickwidth = 2.5
linewidth = 2


def rb():
    pass

def epsb_dot(all_param, sn):
    # xs,ys,zs,rs = sn['PartType5']['Coordinates']
    G = all_param['GravityConstantInternal']
    mb = 1
    mg = sn['PartType0']['Masses']

    # // By convention:
    # // vb = v1-v2
    vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]
    dxb = sn['PartType5']['Coordinates'][0][0] - sn['PartType5']['Coordinates'][1][0]
    dyb = sn['PartType5']['Coordinates'][0][1] - sn['PartType5']['Coordinates'][1][1]
    rb = [dxb, dyb]
    rb_mag = np.sqrt(dxb**2 + dyb**2)
    epsb = 0.5 * np.sqrt(vb[0]**2 + vb[1]**2) - G*mb/rb_mag

    dx1 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][0][0]
    dy1 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][0][1]
    dr1 = np.sqrt(dx1**2 + dy1**2)
    fgrav1 = -1 * G * sn['PartType0']['Masses'] * np.array([dx1, dy1]) * 1./(dr1**3)

    dx2 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][1][0]
    dy2 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][1][1]
    dr2 = np.sqrt(dx2**2 + dy2**2)
    fgrav2 = -1 * G * sn['PartType0']['Masses'] * np.array([dx2, dy2]) * 1./(dr2**3)

    fgrav = fgrav1 + fgrav2

    depsb = vb[0]*fgrav[0]+vb[1]*fgrav[1]

    return(depsb, epsb, fgrav)

def lb_dot(all_param, sn, fgrav):
    vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]
    dxb = sn['PartType5']['Coordinates'][0][0] - sn['PartType5']['Coordinates'][1][0]
    dyb = sn['PartType5']['Coordinates'][0][1] - sn['PartType5']['Coordinates'][1][1]
    rb = [dxb, dyb]
    dlb = rb[0]*fgrav[1] - rb[1]*fgrav[0]
    lb = rb[0]*vb[1] - rb[1]*vb[0]
    return(dlb, lb)


def eb_dot(all_param, sn):
    depsb, epsb, fgrav = epsb_dot(all_param, sn)
    dlb, lb = lb_dot(all_param, sn, fgrav)

    eb = all_param['Eccentricity']

    mbdot = 0
    mb = 1

    if epsb > 0:
        return((1 - eb**2)/(2.*eb) * \
            (2. * mbdot/mb - depsb/epsb - (2. * dlb/lb)))
    else:
        return((1 - eb**2) * 
            (2. * mbdot/mb - depsb/epsb - (2. * dlb/lb)))

def plot(all_fps, figpath, titles, box=20, BoxSize=300, \
    dim='2D', SemiMajorAxis=1, fname='ebdot', contour=False, ratio = False, **kwargs):
    center=[BoxSize/2., BoxSize/2.]

    all_param_fp = open(all_fps[0] + '/../' + 'all_param.pkl', 'rb')
    all_param = pkl.load(all_param_fp)
    all_param['BoxSize'] = BoxSize

    i_initial, i_final = misc.snap_inds(all_fps, kwargs)

    # if 'vmax' not in kwargs and 'vmin' not in kwargs:
    #     [vmin, vmax] = rad.find_ybounds_multiple(all_fps, 'Viscosity', i_initial, i_final)
    # else:
    vmax = kwargs['vmax']
    vmin = kwargs['vmin']

    for i in range(i_initial, i_final+1):
        ''' Make figure and axes '''
        n_ax = len(all_fps)
        fig = plt.figure(figsize = (n_ax*figwidth, figheight))
        sn_time = gadget.Simulation(all_fps[0] + '/snap_%03d.hdf5' %i)
        # fig.suptitle("t = %.2f" %(sn_time.Time.value))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        all_ebdot = {}

        for k,fp in enumerate(all_fps):
            ''' First load snapshots '''
            sn = h5.File(fp +  '/snap_%03d.hdf5' %i, 'r')
            sn_sim = gadget.Simulation(fp + '/snap_%03d.hdf5' %i)
            ax = fig.add_subplot(1, n_ax, k+1)
            ax.set_title(titles[k] + r'$, t = %.2f\, P_{\rm b}$' %(sn_time.Time.value/(2*np.pi)), fontsize=fs)
            ''' Get x and y coords of snapshot '''
            x,y,z,r = misc.x_coords(sn, dim=dim, BoxSize=BoxSize)

            ''' Calculate ebdot for each gas cell in domain '''
            all_ebdot['ebdot_%d' %k] = eb_dot(all_param, sn)
            #vmin = vmin, vmax = vmax, 
            im = sn_sim.plot_Aslice(all_ebdot['ebdot_%d' %k], center = center, \
                                    contour=contour, colorbar = False, axes = ax, box = box, \
                                    norm=colors.SymLogNorm(linthresh=1.e-10, linscale=0.03, vmin=-1.e-6, vmax=1.e-6), \
                                    vmin=-1.e-6, vmax=1.e-6, \
                                    cmap='RdBu_r')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)

            ''' Add a marker for the sink particles '''
            for sink_no in np.arange(0,len(sn['PartType5']['ParticleIDs'])):
                marker_coords = sn['PartType5']['Coordinates'][sink_no][0:2]
                print("marker_coords = ", marker_coords)
                circle = plt.Circle((marker_coords[0], marker_coords[1]), 0.05, linewidth=3, color='r', fill=False)
                ax.add_patch(circle)
                #ax.Circle((marker_coords[0], marker_coords[1]), 0.05, linewidth=3, color='r', fill=False)
                #ax.plot(marker_coords[0], marker_coords[1], 'o', markersize = 0.5, color='r')
            plt.tight_layout()

        misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
        
        plt.savefig(figpath + fname + '_%03d.png' %i)
        plt.close()

    # misc.make_video(figpath, fname)