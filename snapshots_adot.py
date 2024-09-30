import h5py as h5
import glob
import gadget
from gadget import calcGrid
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
# setting path
sys.path.append('../')
import misc as misc
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
import misc
import torques as tq
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import inspector_gadget_sliceplots as igs
import eggleton_roche as er

#SET A COLORBLIND COMPATIBLE SCHEME FOR ENTIRE PAPER
import seaborn as sns
palette = sns.color_palette('colorblind')
palette2 = sns.color_palette('bright')
all_colors = {
'e=0.00': palette[0],
'e=0.10': palette[9],
'e=0.20': palette[1],
'e=0.30': palette[7],
'e=0.40': palette[2],
'e=0.60': palette[3],
'e=0.80': palette[4],
}

rs_colors = {
    'rs=0.030': 'red',
    'rs=0.010': 'orange',
    'rs=0.005': 'green',
}

figwidth = 16
figheight = 14
lw = 15
ticksize = 40
tickwidth = 6

def plot_figure(all_fps, labels, figpath, torque_contributions, \
                all_accrad, snapshot_ids, \
                vmin = 1.e-7, vmax = 1.e-4, contour=False, \
                fname='torques_snapshot_sync_rs_param_study', \
                primary = True, secondary = False, \
                n_frames = 10, widths=[1,1,1], **kwargs):

    fs = 20 * len(snapshot_ids)
    all_be_arrays = []
    all_param_dicts = []

    for fp in all_fps:
        tmin = gadget.Simulation(fp +  '/snap_%03d.hdf5' %snapshot_ids[0]).Time.value
        tmax = gadget.Simulation(fp +  '/snap_%03d.hdf5' %snapshot_ids[-1]).Time.value
        print("line 89 in snapshots_adot")
        binary_evolution_array = be.load_binary_evolution(fp, tmin=tmin, tmax=tmax)
        print("line 91 in snapshots_adot")
        all_param = misc.load_param_txt(fp)
        print("line 93 in snapshots_adot")
        all_be_arrays.append(binary_evolution_array)
        all_param_dicts.append(all_param)
        print("line 96 in snapshots_adot")

    t_orb = misc.T_orbit(a = all_param['SemiMajorAxis'])
    ax_acc_ylim_lower = []
    ax_acc_ylim_upper = []
    ax_abdot_ylim_lower = []
    ax_abdot_ylim_upper = []

    fig = plt.figure(figsize=(len(snapshot_ids)*figwidth, (len(all_accrad)+len(torque_contributions))*figheight))
    gs = gridspec.GridSpec((len(all_accrad)+len(all_accrad)), len(snapshot_ids),\
                            height_ratios = [1]*(len(all_accrad)+len(all_accrad)))

    all_ax_abdot = []
    for k,accrad in enumerate(all_accrad):
        if k > 0:
            all_ax_abdot.append(plt.subplot(gs[-(k+1), :], sharex = all_ax_abdot[0]))
        else:
            all_ax_abdot.append(plt.subplot(gs[-(k+1), :]))
    #all_ax_abdot.reverse()
    print("all_ax_abdot = ", all_ax_abdot)

    for i,(fp,binary_evolution_array,all_param,label,accrad) in enumerate(zip(all_fps,all_be_arrays,all_param_dicts,labels,all_accrad)):
        for j,snap_id in enumerate(snapshot_ids):
            ax = plt.subplot(gs[i, j])
            sn = gadget.Simulation(fp + '/', snap_id)
            kwargs['i_initial'] = snap_id
            kwargs['i_final'] = snap_id
            kwargs['all_param'] = all_param_dicts[i]
            all_param['AccretionRadius'] = accrad
            sigma_0 = all_param['rho_0']
            box = widths[i]
            label = labels[i]

            if snap_id == snapshot_ids[-1]:
                ax, im_p = igs.image_snaps(fp, figpath, ax = ax, title = '', \
                          BoxSize=all_param['BoxSize'], \
                          sink_particles = True, save= False, video_only = False, video = False,\
                          n_frames=n_frames, vmin = vmin/sigma_0, vmax = vmax/sigma_0, contour = contour, \
                          colorbar_fraction = 0.05, \
                          SemiMajorAxis=all_param['SemiMajorAxis'], \
                          centered_on_primary=primary, \
                          centered_on_secondary=secondary, \
                          show_colorbar = False, \
                          fontsize=fs, sigma_0 = sigma_0, show_time=True, \
                          box = box, return_image = True, **kwargs)
            else:
                ax = igs.image_snaps(fp, figpath, ax = ax, title = '', \
                          BoxSize=all_param['BoxSize'], \
                          sink_particles = True, save= False, video_only = False, video = False,\
                          n_frames=n_frames, vmin = vmin/sigma_0, vmax = vmax/sigma_0, contour = contour, \
                          colorbar_fraction = 0.05, \
                          SemiMajorAxis=all_param['SemiMajorAxis'], \
                          centered_on_primary=primary, \
                          centered_on_secondary=secondary, \
                          show_colorbar = False, \
                          fontsize=fs, sigma_0 = sigma_0, show_time=True, \
                          box = box, return_image = False, **kwargs)
                    
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xlabel('')
            ax.set_ylabel('')
            # plt.setp(ax.get_xticklabels(which='both'), fontsize=0.7*fs)
            # plt.setp(ax.get_yticklabels(which='both'), fontsize=0.7*fs)
            ax.tick_params(axis='both', which='both', direction='inout', size=ticksize, width=tickwidth)
            plt.tight_layout()

            if snap_id == snapshot_ids[0]:
                ax.set_ylabel(r'$r_{\rm s} = $' + '%.3f' %accrad, rotation=90, fontsize=1.25*fs)

            if snap_id == snapshot_ids[-1]:
                """ ADD COLORBAR TO SECOND AXIS """
                cb = plt.colorbar(im_p, pad=0.1, ax=ax, \
                                  fraction=0.10, aspect = 20)
                cb.set_label(label=r'$\Sigma/\Sigma_0$',size=1.2*fs,weight='bold')
                cb.ax.tick_params(labelsize=fs, which='both', size=ticksize, width=tickwidth)
                """---------------PLOT PROJECTIONS FIRST--------------"""
            plt.tight_layout()

        #for accrad,ax_abdot in zip(all_accrad,all_ax_abdot):
        ax_abdot = all_ax_abdot[-(i+1)]
        linestyles = ['-', '--']
        for torque_contribution,ls in zip(torque_contributions,linestyles):
            color = rs_colors['rs=%.3f' %accrad]
            """---------------PLOT TIME-SERIES OF \dot{a} --------------"""
            if torque_contribution == 'sum_grav':
                torque_label = r'$ \rm{g} $'
            if torque_contribution == 'sum_grav_a':
                torque_label = r'$ \rm{g,r>a} $'
            if torque_contribution == 'sum_grav_acc':
                torque_label = r'$ \rm{g+a} $'

            """ FIND INDICES WHERE ARRAY QUANTITIES ARE WITHIN DESIRED RANGE OF tmin<=t=<tmax """
            ind_tmin = binary_evolution_array['t'] >= tmin
            ind_tmax = binary_evolution_array['t'] <= tmax
            ind_t = ind_tmin & ind_tmax
            """ FIND INDICES WHERE ARRAY QUANTITIES ARE WITHIN DESIRED RANGE OF tmin<=t=<tmax """
            ab_dot_ab = binary_evolution_array['ab_dot_ab_%s' %torque_contribution][ind_t]
            eb_dot = binary_evolution_array['eb_dot_%s' %torque_contribution][ind_t]
            t = binary_evolution_array['t'][ind_t]            
            # + ', mean=%.2e' %(np.mean(ab_dot_ab))
            ax_abdot.plot(t/(2*np.pi), ab_dot_ab, label=torque_label, \
                          linewidth=lw, color=color, linestyle = ls, alpha = 0.6)
            ax_abdot.set_ylabel(r'$\dot{a}_{\rm b}/a_{\rm b}$', fontsize=fs)
            ax_abdot.grid(which='both', color='k', linestyle='-', linewidth=0.5)
            x_loc_a2 = 0.9
            y_loc_a2 = 0.2
            title_out = r'$r_{\rm s} = $' + '%.3f' %accrad
            t = ax_abdot.text(x_loc_a2, y_loc_a2, title_out,\
            horizontalalignment='center', \
            verticalalignment='center', \
            transform = ax_abdot.transAxes, fontsize=fs)
            t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))
            ax_abdot.tick_params(axis='both', which='both', direction='inout', size=ticksize, width=tickwidth)
            plt.tight_layout()

        for ax in all_ax_abdot:
            ax.legend(fontsize=fs, ncol=1, loc='lower left')
            plt.setp(ax.get_yticklabels(which='both'), fontsize=fs)
            if ax == all_ax_abdot[0]:
                ax.set_xlim([tmin/t_orb, tmax/t_orb])
                plt.setp(ax.get_xticklabels(which='both'), fontsize=fs)
                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                ax.set_xlabel(r'$t/P_{\rm b}$', fontsize=fs)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)
            plt.tight_layout()

    plt.tight_layout()
    fig.savefig(figpath + fname + '.png')
    plt.close()


if __name__ == "__main__":
    alpha = 0.1
    accrad = 0.03
    sinkrate = 0.50
    phi = 0.00
    e = 0.8
    q = 1.0
    h = 0.10
    dim = 2
    res = 5.e-10
    box = 1 
    snapshot_ids = 12000 + np.array([20, 45, 50, 55, 80]) #np.array([45, 48, 50, 52, 55])
    primary = False
    secondary = True
    #[800, 825, 850, 875, 900]

    fp_root = '/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_eccentric_accretion_paper_resub_videos/'
    figpath = '/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/figures/torque_paper/binary_evolution_csd_mass_mdot_snapshot_video/'
    if not os.path.exists(figpath):
        print("Path doesn't exist, making new filepath: ", figpath)
        os.makedirs(figpath)

    fname = 'snapshots_adot_eb%.2f_qb%.2f' %(e,q)

    if primary:
        fname += '_primary'
    if secondary:
        fname += '_secondary'

    all_fps = []
    labels = []
    all_accrad = [0.03, 0.005]
    for accrad in all_accrad:
        fp = fp_root + 'res=%.3e/alpha=%.3f/dim=%d/e=%.3f/q=%.3f/aspect_ratio=%.3f/phi=%.3f/AccretionRadius=%.3f/SinkRate=%.3f/output/' \
                                %(res,alpha,dim,e,q,h,phi,accrad,sinkrate)
        all_fps.append(fp)
        labels.append(r'$r_{\rm s} = $' + '%.3f' %accrad)
        #fname += '_%.3f' %accrad
    
    torque_contributions = ['sum_grav_acc', 'sum_grav']

    widths = [0.5]*len(all_accrad)

    # for torque_contribution in torque_contributions:
    #     fname += '_%s' %torque_contribution
    
    plot_figure(all_fps, labels, figpath, torque_contributions, \
                all_accrad, snapshot_ids, \
                vmin = 1.e-7, vmax = 1.e-4, contour=False, \
                fname=fname, primary=primary, secondary=secondary, \
                n_frames = 10, widths=widths)

