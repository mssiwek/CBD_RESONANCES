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
import accretion_theory as act


def calc_mbdot(all_param):
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
    epsb = - G*mb/(2*all_param['SemiMajorAxis'])

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

    mbdot_norm = calc_mbdot(all_param)
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

    mbdot = calc_mbdot(all_param)
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
    # mbdot = calc_mbdot(all_param)
    mbdot = 0
    mb = all_param['CentralMass']
    depsb, epsb = calc_epsbdot(all_param, sn)
    ab = all_param['SemiMajorAxis']

    abdot = mbdot/mb - depsb/epsb 

    return(abdot, ab)

def calc_sum_abdot(all_param, sn):
    depsb, epsb = calc_epsbdot(all_param, sn)

    eb = all_param['Eccentricity']
    rb_mag = calc_rb_mag(sn)

    mbdot = calc_mbdot(all_param)
    mb = all_param['CentralMass']

    #sum all the depsb and dlb
    depsb_sum = sum(depsb) - mbdot/rb_mag

    abdot = mbdot/mb - depsb_sum/epsb 

    return(abdot)