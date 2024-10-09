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

figwidth = 7
figheight = 6
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2

def calc_mbdot(all_param, sn):
    acc = all_param['acc_txt']
    snap_time = all_param['snap_time']
    acc_t_abs = np.abs(acc['t']-snap_time)
    #find the index of the timestamp in accretion.txt closest to the snapshot time
    t_idx = acc_t_abs.loc[(acc_t_abs == acc_t_abs.min())].index[0]
    return(acc['mdot'][t_idx])

def calc_rb_mag(sn):
    dxb = sn['PartType5']['Coordinates'][0][0] - sn['PartType5']['Coordinates'][1][0]
    dyb = sn['PartType5']['Coordinates'][0][1] - sn['PartType5']['Coordinates'][1][1]
    rb_mag = np.sqrt(dxb**2 + dyb**2)
    return(rb_mag)


def calc_fgrav(all_param, sn):
    G = all_param['GravityConstantInternal']
    mb = 1

    # By convention:
    # vb = v1-v2
    # fgrav = fgrav1 - fgrav2

    dx1 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][0][0]
    dy1 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][0][1]
    dr1 = np.sqrt(dx1**2 + dy1**2)
    fgrav1 = -1 * G * sn['PartType0']['Masses'] * np.array([dx1, dy1]) * 1./(dr1**3)

    dx2 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][1][0]
    dy2 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][1][1]
    dr2 = np.sqrt(dx2**2 + dy2**2)
    fgrav2 = -1 * G * sn['PartType0']['Masses'] * np.array([dx2, dy2]) * 1./(dr2**3)

    fgrav = fgrav1 - fgrav2

    return(fgrav, fgrav1, fgrav2)


def calc_epsbdot(all_param, sn):
    #not including the term (- G * Mbdot)/rb
    G = all_param['GravityConstantInternal']
    mb = 1
    # xs,ys,zs,rs = sn['PartType5']['Coordinates']
    fgrav, *_ = calc_fgrav(all_param, sn)

    vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]

    rb_mag = calc_rb_mag(sn)
    # epsb = 0.5 * np.sqrt(vb[0]**2 + vb[1]**2) - G*mb/rb_mag
    epsb = - G*mb/(2*rb_mag)

    depsb = vb[0]*fgrav[0]+vb[1]*fgrav[1] 

    return(depsb, epsb)

def calc_lbdot(all_param, sn):
    fgrav, *_ = calc_fgrav(all_param, sn)

    vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]
    dxb = sn['PartType5']['Coordinates'][0][0] - sn['PartType5']['Coordinates'][1][0]
    dyb = sn['PartType5']['Coordinates'][0][1] - sn['PartType5']['Coordinates'][1][1]
    rb = [dxb, dyb]
    dlb = rb[0]*fgrav[1] - rb[1]*fgrav[0]
    lb = np.sqrt(all_param['GravityConstantInternal'] * all_param['CentralMass'] \
                * all_param['SemiMajorAxis'] * (1. - all_param['Eccentricity']**2))
    return(dlb, lb)

def calc_ebdot(all_param, sn):
    depsb, epsb = calc_epsbdot(all_param, sn)
    dlb, lb = calc_lbdot(all_param, sn)

    eb = all_param['Eccentricity']

    # Here we are calculating ebdot _for each cell_.
    # It only makes sense to add mbdot to get net ebdot 
    # after summing ebdot over entire domain.

    mbdot_norm = calc_mbdot(all_param, sn)
    mbdot = 0
    mb = all_param['CentralMass']

    if eb > 0:
        ebdot = (((1 - eb**2)/(2.*eb)) * \
            (2. * mbdot/mb - depsb/epsb - (2. * dlb/lb))) * mb/mbdot_norm
    else:
        ebdot = ((1 - eb**2) * \
            (2. * mbdot/mb - depsb/epsb - (2. * dlb/lb))) * mb/mbdot_norm
    
    return(ebdot, eb)

def calc_sum_ebdot(all_param, sn):
    depsb, epsb = calc_epsbdot(all_param, sn)
    dlb, lb = calc_lbdot(all_param, sn)
    dlb_sum = sum(dlb)

    eb = all_param['Eccentricity']
    rb_mag = calc_rb_mag(sn)

    mbdot = calc_mbdot(all_param, sn)
    mb = all_param['CentralMass']

    #sum all the depsb and dlb
    depsb_sum = sum(depsb) - mbdot/rb_mag
    dlb_sum = sum(dlb)

    if eb > 0:
        ebdot = (((1 - eb**2)/(2.*eb)) * \
            (2. * mbdot/mb - depsb_sum/epsb - (2. * dlb_sum/lb))) * mb/mbdot
    else:
        ebdot = ((1 - eb**2) * \
            (2. * mbdot/mb - depsb_sum/epsb - (2. * dlb_sum/lb))) * mb/mbdot
    
    return(ebdot)



def calc_abdot(all_param, sn):
    # mbdot = calc_mbdot(all_param, sn)
    mbdot = 0
    mb = all_param['CentralMass']
    depsb, epsb = calc_epsbdot(all_param, sn)
    ab = all_param['SemiMajorAxis']

    abdot = mbdot/mb - depsb/epsb 

    return(abdot, ab)


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
        n_ax = len(all_fps)
        fig = plt.figure(figsize = (n_ax*figwidth, figheight))
        sn_time = gadget.Simulation(all_fps[0] + '/snap_%03d.hdf5' %i)
        all_param['snap_time'] = sn_time.Time.value
        # fig.suptitle("t = %.2f" %(sn_time.Time.value))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        vals = {}

        for k,(fp,all_param) in enumerate(zip(all_fps,all_params)):
            ''' First load snapshots '''
            sn = h5.File(fp +  '/snap_%03d.hdf5' %i, 'r')
            sn_sim = gadget.Simulation(fp + '/snap_%03d.hdf5' %i)
            ax = fig.add_subplot(1, n_ax, k+1)
            ''' Get x and y coords of snapshot '''
            x,y,z,r = misc.x_coords(sn, dim=dim, BoxSize=BoxSize)

            ''' Calculate param for each gas cell in domain '''
            vals['%s_%d' %(param,k)], *_ = eval('calc_%s(all_param, sn)' %param)
            ''' Calculate the sum of ebdot and abdot '''
            vals['%s_sum_%d' %(param,k)] = eval('calc_sum_%s(all_param, sn)' %param)

            ax.set_title(titles[k] + r'$, t = %.2f\, P_{\rm b}$' %(sn_time.Time.value/(2*np.pi)) + "\n" + \
                        r'$\dot{%s}_{\rm b} = %.2f \dot{M}_{\rm b}/M_{\rm b}$' %(param[0], vals['%s_sum_%d' %(param,k)]), \
                        fontsize=fs)

            v_range = 1.e-2
            cbar_norm_log = colors.SymLogNorm(linthresh=1.e-5, linscale=0.03, \
                                    vmin=-v_range, vmax=v_range)
            cbar_norm_lin = mpl.colors.Normalize(vmin=-v_range, vmax=v_range)
            
            im = sn_sim.plot_Aslice(vals['%s_%d' %(param,k)], center = center, \
                                    contour=contour, colorbar = False, axes = ax, box = box, \
                                    norm=cbar_norm_lin, \
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

        misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
        
        plt.savefig(figpath + fname + '_%03d.png' %i)
        plt.close()

    # misc.make_video(figpath, fname)