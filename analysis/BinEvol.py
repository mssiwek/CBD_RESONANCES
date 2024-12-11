""" 
    THIS CLASS CONSTRUCTS AN OBJECT WHICH CONTAINS 
    BINARY EVOLUTION RATES \dot{a}_{\rm b} and \dot{e}_{\rm b} 
    FOR A GIVEN SET OF PARAMETERS.
    THE CLASS TAKES ON-THE-FLY ANGULAR MOMENTUM AND ENERGY RATES 
    AVAILABLE IN torques_eb_lb.txt FROM SIWEK+23B,
    BUT USER-DEFINED INPUT VALUES CAN ALSO BE SUPPLIED

    THE GOAL IS TO END UP WITH A TIME-SERIES OF;
    LBDOT, EPSBDOT, EBDOT, ABDOT

    created by Magdalena Siwek, last modified 10/2024
"""

import numpy as np  
import misc
import sys
import time
import h5py as h5
import pandas as pd

class BinEvol:
    def __init__(self, SimInit, from_txt=True, from_snap=False, sn=None):
        """ IMPORT THE INFO NEEDED TO COMPUTE ANGULAR MOMENTUM 
        AND ENERGY RATE OF CHANGE (IN CASE OF SNAPSHOT),
        OR GET DLB AND DEPSB DIRECTLY (IF USING torques_eb_lb.txt) """

        if from_txt:
            SimInit = SimInit.evol_from_txt()
        elif from_snap:
            SimInit = SimInit.evol_from_snap(sn)
        else:
            sys.exit("Need to calculate evolution rates either from .txt or snapshot.")
        
        self.SimInit = SimInit
        self.from_txt = from_txt
        self.from_snap = from_snap

        """ LOAD COMMONLY USED CONSTANTS """
        self.G = self.SimInit.all_param['GravityConstantInternal']
        self.mb = self.SimInit.all_param['CentralMass']
        self.ab = self.SimInit.all_param['SemiMajorAxis']
        self.eb = self.SimInit.all_param['Eccentricity']
        self.qb = self.SimInit.all_param['MassRatio']


    def ebdot(self):
        if self.eb > 0:
            ebdot = ((1 - self.eb**2)/(2.*self.eb)) \
                    * (2. * self.SimInit.mbdot/self.mb \
                    - self.SimInit.epsbdot/self.SimInit.epsb \
                    - (2. * self.SimInit.lbdot/self.SimInit.lb))
        else:
            ebdot = (1 - self.eb**2) \
                    * (2. * self.SimInit.mbdot/self.mb \
                    - self.SimInit.epsbdot/self.SimInit.epsb \
                    - (2. * self.SimInit.lbdot/self.SimInit.lb))
        return(ebdot)

    def abdot(self):
        abdot = self.SimInit.mbdot/self.mb \
                - self.SimInit.epsbdot/self.SimInit.epsb
        return(abdot) 
    
    def mbdot(self):
        return(self.SimInit.mbdot)
    
    def mbdot_norm(self):
        return(self.SimInit.mbdot_norm)



""" INITIALIZE THE SIMULATION. RETURNS THE LOCATIONS AND MASSES 
    OF GAS PARTICLES WE WANT TO USE TO COMPUTE ebdot AND abdot 
    IF USING A SNAPSHOT. (snap = True)

    IF USING A TXT FILE, SEARCHES FOR torques_eb_lb.txt TO LOAD
    DLB AND DEPSB.

    RG: LOCATION OF GAS PARTICLE. FORMAT:
    RG[0][X,Y,Z]: PRIMARY
    RG[1][X,Y,Z]: SECONDARY
    
    MB[0]: MASS OF PRIMARY, MB[1]: MASS OF SECONDARY"""
class SimInit:
    def __init__(self, ext, acc=False, f_acc=False, fg_cav=True, **kwargs):
        self.acc = acc
        self.f_acc = f_acc
        self.fg_cav = fg_cav
        self.ext = ext
        start_time = time.time()
        all_param = misc.load_param_txt(self.ext + '/../')
        print("loading all_param.txt took %s seconds ---" % (time.time() - start_time))

        self.all_param = all_param

        self.epsb = - self.all_param['GravityConstantInternal']\
                    * self.all_param['CentralMass']/(2*self.all_param['SemiMajorAxis'])
        self.lb = np.sqrt(self.all_param['GravityConstantInternal']\
                          * self.all_param['CentralMass']\
                          * self.all_param['SemiMajorAxis']\
                          * (1. - self.all_param['Eccentricity']**2))
        
        if self.all_param['InclinationAngle'] > 0:
            print("Inclination Angle = %.2f. Assuming the binary is retrograde, \n and therefore lb is multiplied by -1." %(self.all_param['InclinationAngle']))
            self.lb *= (-1)

    
    def evol_from_txt(self):
        """ THIS MEANS WE ARE SUPPLYING A FILEPATH VIA ext,
            AND GET dlb AND depsb FROM torques_eb_lb FILE """
        
        """ LOAD SIMULATION OUTPUT AND PARAMETERS """
        """ ------------------------------------- """
        """ ------------------------------------- """
        start_time = time.time()
        tqfl = self.read_torques_txt(orbit_averaged=True, n_orbit_average = 1, \
                                     raw=False, dot_eb = True, old_sims=True)
        print("loading torques_eb_lb.txt took %s seconds ---" % (time.time() - start_time))

        """ ------------------------------------- """
        """ ------------------------------------- """
        """ LOAD SIMULATION OUTPUT AND PARAMETERS """

        start_time = time.time()
        """ LOADING BINARY POSITION AND VELOCITY """
        """ ------------------------------------ """
        """ ------------------------------------ """
        self.rb_mag = tqfl['r_b_mag']
        self.rb = tqfl['r_b']
        self.vb = tqfl['v_b']
        self.t = tqfl['t']
        """ ------------------------------------ """
        """ ------------------------------------ """
        """ LOADING BINARY POSITION AND VELOCITY """



        """ ------------------------------------------------------ """
        """ ------------------------------------------------------ """
        """ GETTING ANGULAR MOMENTUM INCREMENTAL CHANGE (delta lb) """
        if self.fg_cav:
            self.dlb1 = tqfl['dlb_sum_grav_1_2']
            self.dlb2 = tqfl['dlb_sum_grav_2_2']
        else:
            """ IF fg_cav == False, ONLY GAS CELLS 
            WITH r>a_b ARE INCLUDED IN CALCULATION:
                f_ext = f_grav[r>a_b] """
            self.dlb1 = tqfl['dlb_sum_grav_a_1_2']
            self.dlb2 = tqfl['dlb_sum_grav_a_2_2']
        if self.f_acc:
            if not self.fg_cav:
                exit("Cannot include accretion torques if we choose to exclude r<a_b gas cells!")
            else:
                """ f_acc REFERS TO THE ACCRETION OF LINEAR MOMENTUM, I.E. 
                    IF f_acc == True, WE ADD THIS TO THE 
                    EXTERNAL FORCE APPLIED TO BINARY:
                    f_ext = f_grav + f_acc """
                self.dlb1 += tqfl['dlb_sum_acc_1_2']
                self.dlb2 += tqfl['dlb_sum_acc_2_2']
        """ GETTING ANGULAR MOMENTUM INCREMENTAL CHANGE (delta lb) """
        """ ------------------------------------------------------ """
        """ ------------------------------------------------------ """



        """ --------------------------------- """
        """ --------------------------------- """
        """ GETTING ENERGY INCREMENTAL CHANGE """
        if not self.fg_cav:
            """ EXCLUDING cavity """
            """ f_ext = f_grav ONLY """
            self.depsb1 = tqfl['depsb_sum_grav_a_no_mdot_1']
            self.depsb2 = tqfl['depsb_sum_grav_a_no_mdot_2']
        else:
            """ INCLUDING cavity """
            if self.f_acc:
                """ f_acc REFERS TO THE ACCRETION OF LINEAR MOMENTUM, I.E. 
                    IF f_acc == True, WE ADD THIS TO THE 
                    EXTERNAL FORCE APPLIED TO BINARY:
                    f_ext = f_grav + f_acc """
                self.depsb1 = tqfl['depsb_sum_grav_acc_no_mdot_1']
                self.depsb2 = tqfl['depsb_sum_grav_acc_no_mdot_2']
            else:
                """ f_ext = f_grav ONLY """
                self.depsb1 = tqfl['depsb_sum_grav_no_mdot_1']
                self.depsb2 = tqfl['depsb_sum_grav_no_mdot_2']
        """ GETTING ENERGY INCREMENTAL CHANGE """
        """ --------------------------------- """
        """ --------------------------------- """


        if self.all_param['MassRatio'] == 1:
            """ FOR EQUAL MASS RATIOS, THE ORDER DEPENDS ON THE ID VALUE """
            if tqfl['id_1'][0] > tqfl['id_2'][0]:
                """ IN THIS CASE, "PRIMARY" IS DEFINED AS id_1
                    JUST NEED TO BE CONSISTENT WITH THIS THROUGHOUT """
                s1 = +1
                s2 = -1
            else:
                """ IN THIS CASE, "PRIMARY" IS DEFINED AS id_2 
                    JUST NEED TO BE CONSISTENT WITH THIS THROUGHOUT """
                s1 = -1
                s2 = +1
        else:
            s1 = +1
            s2 = -1
        
        """ CALCULATE TOTAL INCREMENTAL CHANGE 
            IN ANGULAR MOMENTUM AND ENERGY """
        """ --------------------------------- """
        """ --------------------------------- """
        self.dlb = s1 * self.dlb1 + s2 * self.dlb2
        self.depsb = s1 * self.depsb1 + s2 * self.depsb2
        """ --------------------------------- """
        """ --------------------------------- """
        """ CALCULATE TOTAL INCREMENTAL CHANGE 
            IN ANGULAR MOMENTUM AND ENERGY """
        
        if not self.fg_cav:
            if self.acc:
                """ NOW WE ARE ALSO ADDING EFFECT OF ACCRETION TO THE 
                    RATE OF CHANGE OF EPSILON (BINARY SPECIFIC ENERGY) """
                self.depsb += - tqfl['dm']/tqfl['r_b_mag']
                self.dmb = tqfl['dm']
            else:
                self.dmb = 0
                self.mbdot = 0
        else:
            self.dmb = 0
            self.mbdot = 0

        self.epsbdot = np.zeros(np.shape(tqfl['t']))
        self.epsbdot[1:] = self.depsb[1:]/np.diff(tqfl['t'])

        self.lbdot = np.zeros(np.shape(tqfl['t']))
        self.lbdot[1:] = self.dlb[1:]/np.diff(tqfl['t'])

        self.mbdot_norm = np.zeros(np.shape(tqfl['t']))
        self.mbdot_norm[1:] = tqfl['dm'][1:]/np.diff(tqfl['t'])
        self.mbdot_norm[0] = self.mbdot_norm[1]
        
        """ dmb to normalize! """
        self.dmb_norm = tqfl['dm']
        print("remainder of evol_from_txt took %s seconds ---" % (time.time() - start_time))

        return(self)
        

    def evol_from_snap(self, sn):
        """Note: ebdot and abdot maps do not include accretion contributions."""
        """ MUST GET LOCATIONS AND MASSES 
            OF GAS DISTRIBUTION WE WANT TO USE 
            TO CALCULATE BINARY EVOLUTION """
        
        if sn == None:
            exit("Must supply snapshot to calculate evolution from snapshot!")
        
        """ Gravitational acceleration on binary from gas """
        G = self.all_param['GravityConstantInternal']
        # By convention:
        # vb = v1-v2
        # fgrav = fgrav1 - fgrav2

        """ LOADING BINARY AND GAS COORDINATES """
        """ ---------------------------------- """
        """ ---------------------------------- """
        dx1 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][0][0]
        dy1 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][0][1]
        dr1 = np.sqrt(dx1**2 + dy1**2)
        dx2 = sn['PartType0']['Coordinates'][:,0] - sn['PartType5']['Coordinates'][1][0]
        dy2 = sn['PartType0']['Coordinates'][:,1] - sn['PartType5']['Coordinates'][1][1]
        dr2 = np.sqrt(dx2**2 + dy2**2)
        """ ---------------------------------- """
        """ ---------------------------------- """
        """ LOADING BINARY AND GAS COORDINATES """
        

        """ GRAVITATIONAL ACCELERATION GAS-BINARY """
        """ ------------------------------------- """
        """ ------------------------------------- """
        fgrav1 = -1 * G * sn['PartType0']['Masses'] * np.array([dx1, dy1]) * 1./(dr1**3)
        fgrav2 = -1 * G * sn['PartType0']['Masses'] * np.array([dx2, dy2]) * 1./(dr2**3)
        self.fgrav = fgrav1 - fgrav2
        """ ------------------------------------- """
        """ ------------------------------------- """
        """ GRAVITATIONAL ACCELERATION GAS-BINARY """


        """ --------------------------------------- """
        """ --------------------------------------- """
        """ GETTING ANGULAR MOMENTUM RATE OF CHANGE """
        dxb = sn['PartType5']['Coordinates'][0][0] - sn['PartType5']['Coordinates'][1][0]
        dyb = sn['PartType5']['Coordinates'][0][1] - sn['PartType5']['Coordinates'][1][1]
        rb = [dxb, dyb]
        self.lbdot = rb[0]*self.fgrav[1] - rb[1]*self.fgrav[0]
        """ GETTING ANGULAR MOMENTUM RATE OF CHANGE """
        """ --------------------------------------- """
        """ --------------------------------------- """
    

        """ ----------------------------- """
        """ ----------------------------- """
        """ GETTING ENERGY RATE OF CHANGE """
        vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]
        self.epsbdot = vb[0]*self.fgrav[0]+vb[1]*self.fgrav[1]
        """ GETTING ENERGY RATE OF CHANGE """
        """ ----------------------------- """
        """ ----------------------------- """

        self.dmb = 0
        self.mbdot = 0
        
        return(self)
    

    def read_torques_txt(self, orbit_averaged=True, n_orbit_average = 1, \
                         raw=False, dot_eb = True, old_sims=True, **kwargs):

        fp = self.ext

        if 'fname' in kwargs:
            fname = kwargs['fname']
        else:
            fname = 'torques_eb_lb'
        all_param = misc.load_param_txt(fp)

        if (all_param['Eccentricity'] == 0.200 or all_param['Eccentricity'] == 0.400) and old_sims:
            columns = ['t', \
                    'id_1', 'id_2', \
                    'mass_1', 'mass_2', \
                    'sink_1_pos_0', 'sink_1_pos_1', 'sink_1_pos_2', \
                    'sink_2_pos_0', 'sink_2_pos_1', 'sink_2_pos_2', \
                    'sink_1_vel_0', 'sink_1_vel_1', 'sink_1_vel_2', \
                    'sink_2_vel_0', 'sink_2_vel_1', 'sink_2_vel_2', \
                    'dm_1', 'dm_2', \
                    'dp_1_0', 'dp_1_1', 'dp_1_2', \
                    'dp_2_0', 'dp_2_1', 'dp_2_2', \
                    'dspin_1_0', 'dspin_1_1', 'dspin_1_2',
                    'dspin_2_0', 'dspin_2_1', 'dspin_2_2',
                    'f_inst_grav_1_0', 'f_inst_grav_1_1', 'f_inst_grav_1_2', \
                    'f_inst_grav_2_0', 'f_inst_grav_2_1', 'f_inst_grav_2_2', \
                    'f_inst_grav_a_1_0', 'f_inst_grav_a_1_1', 'f_inst_grav_a_1_2', \
                    'f_inst_grav_a_2_0', 'f_inst_grav_a_2_1', 'f_inst_grav_a_2_2', \
                    'f_inst_grav_acc_1_0', 'f_inst_grav_acc_1_1', 'f_inst_grav_acc_1_2', \
                    'f_inst_grav_acc_2_0', 'f_inst_grav_acc_2_1', 'f_inst_grav_acc_2_2', \
                    'dlb_sum_grav_1_0', 'dlb_sum_grav_1_1', 'dlb_sum_grav_1_2', \
                    'dlb_sum_grav_2_0', 'dlb_sum_grav_2_1', 'dlb_sum_grav_2_2', \
                    'dlb_sum_grav_a_1_0', 'dlb_sum_grav_a_1_1', 'dlb_sum_grav_a_1_2', \
                    'dlb_sum_grav_a_2_0', 'dlb_sum_grav_a_2_1', 'dlb_sum_grav_a_2_2', \
                    'dlb_sum_acc_1_0', 'dlb_sum_acc_1_1', 'dlb_sum_acc_1_2', \
                    'dlb_sum_acc_2_0', 'dlb_sum_acc_2_1', 'dlb_sum_acc_2_2', \
                    'depsb_sum_grav_1', 'depsb_sum_grav_2', \
                    'depsb_sum_grav_no_mdot_1', 'depsb_sum_grav_no_mdot_2', \
                    'depsb_sum_grav_acc_no_mdot_1', 'depsb_sum_grav_acc_no_mdot_2', \
                    'depsb_sum_grav_acc_1', 'depsb_sum_grav_acc_2', \
                    'depsb_sum_grav_a_1', 'depsb_sum_grav_a_2', \
                    'depsb_sum_acc_1', 'depsb_sum_acc_2']

            columns_bh_1 = ['t', \
                            'id_1',\
                            'mass_1',\
                            'sink_1_pos_0', 'sink_1_pos_1', 'sink_1_pos_2', \
                            'sink_1_vel_0', 'sink_1_vel_1', 'sink_1_vel_2', \
                            'dm_1',\
                            'dp_1_0', 'dp_1_1', 'dp_1_2', \
                            'dspin_1_0', 'dspin_1_1', 'dspin_1_2',
                            'f_inst_grav_1_0', 'f_inst_grav_1_1', 'f_inst_grav_1_2', \
                            'f_inst_grav_a_1_0', 'f_inst_grav_a_1_1', 'f_inst_grav_a_1_2', \
                            'f_inst_grav_acc_1_0', 'f_inst_grav_acc_1_1', 'f_inst_grav_acc_1_2', \
                            'dlb_sum_grav_1_0', 'dlb_sum_grav_1_1', 'dlb_sum_grav_1_2', \
                            'dlb_sum_grav_a_1_0', 'dlb_sum_grav_a_1_1', 'dlb_sum_grav_a_1_2', \
                            'dlb_sum_acc_1_0', 'dlb_sum_acc_1_1', 'dlb_sum_acc_1_2', \
                            'depsb_sum_grav_1', \
                            'depsb_sum_grav_no_mdot_1', \
                            'depsb_sum_grav_acc_no_mdot_1',\
                            'depsb_sum_grav_acc_1', \
                            'depsb_sum_grav_a_1',\
                            'depsb_sum_acc_1']
            columns_bh_2 = ['t', \
                            'id_2',\
                            'mass_2',\
                            'sink_2_pos_0', 'sink_2_pos_1', 'sink_2_pos_2', \
                            'sink_2_vel_0', 'sink_2_vel_1', 'sink_2_vel_2', \
                            'dm_2',\
                            'dp_2_0', 'dp_2_1', 'dp_2_2', \
                            'dspin_2_0', 'dspin_2_1', 'dspin_2_2',
                            'f_inst_grav_2_0', 'f_inst_grav_2_1', 'f_inst_grav_2_2', \
                            'f_inst_grav_a_2_0', 'f_inst_grav_a_2_1', 'f_inst_grav_a_2_2', \
                            'f_inst_grav_acc_2_0', 'f_inst_grav_acc_2_1', 'f_inst_grav_acc_2_2', \
                            'dlb_sum_grav_2_0', 'dlb_sum_grav_2_1', 'dlb_sum_grav_2_2', \
                            'dlb_sum_grav_a_2_0', 'dlb_sum_grav_a_2_1', 'dlb_sum_grav_a_2_2', \
                            'dlb_sum_acc_2_0', 'dlb_sum_acc_2_1', 'dlb_sum_acc_2_2', \
                            'depsb_sum_grav_2', \
                            'depsb_sum_grav_no_mdot_2', \
                            'depsb_sum_grav_acc_no_mdot_2',\
                            'depsb_sum_grav_acc_2', \
                            'depsb_sum_grav_a_2',\
                            'depsb_sum_acc_2']
        else:
            columns = ['t', \
                        'id_1', 'id_2', \
                        'mass_1', 'mass_2', \
                        'sink_1_pos_0', 'sink_1_pos_1', 'sink_1_pos_2', \
                        'sink_2_pos_0', 'sink_2_pos_1', 'sink_2_pos_2', \
                        'sink_1_vel_0', 'sink_1_vel_1', 'sink_1_vel_2', \
                        'sink_2_vel_0', 'sink_2_vel_1', 'sink_2_vel_2', \
                        'dm_1', 'dm_2', \
                        'dp_1_0', 'dp_1_1', 'dp_1_2', \
                        'dp_2_0', 'dp_2_1', 'dp_2_2', \
                        'dspin_1_0', 'dspin_1_1', 'dspin_1_2',
                        'dspin_2_0', 'dspin_2_1', 'dspin_2_2',
                        'f_inst_grav_1_0', 'f_inst_grav_1_1', 'f_inst_grav_1_2', \
                        'f_inst_grav_2_0', 'f_inst_grav_2_1', 'f_inst_grav_2_2', \
                        'f_inst_grav_a_1_0', 'f_inst_grav_a_1_1', 'f_inst_grav_a_1_2', \
                        'f_inst_grav_a_2_0', 'f_inst_grav_a_2_1', 'f_inst_grav_a_2_2', \
                        'f_inst_grav_acc_1_0', 'f_inst_grav_acc_1_1', 'f_inst_grav_acc_1_2', \
                        'f_inst_grav_acc_2_0', 'f_inst_grav_acc_2_1', 'f_inst_grav_acc_2_2', \
                        'dlb_sum_grav_1_0', 'dlb_sum_grav_1_1', 'dlb_sum_grav_1_2', \
                        'dlb_grav_sum_2_0', 'dlb_grav_sum_2_1', 'dlb_sum_grav_2_2', \
                        'dlb_sum_grav_a_1_0', 'dlb_sum_grav_a_1_1', 'dlb_sum_grav_a_1_2', \
                        'dlb_sum_grav_a_2_0', 'dlb_sum_grav_a_2_1', 'dlb_sum_grav_a_2_2', \
                        'dlb_sum_acc_1_0', 'dlb_sum_acc_1_1', 'dlb_sum_acc_1_2', \
                        'dlb_sum_acc_2_0', 'dlb_sum_acc_2_1', 'dlb_sum_acc_2_2', \
                        'depsb_sum_grav_no_mdot_1', 'depsb_sum_grav_no_mdot_2', \
                        'depsb_sum_grav_acc_no_mdot_1', 'depsb_sum_grav_acc_no_mdot_2', \
                        'depsb_sum_grav_a_no_mdot_1', 'depsb_sum_grav_a_no_mdot_2', \
                        'depsb_sum_acc_no_mdot_1', 'depsb_sum_acc_no_mdot_2']
            columns_bh_1 = ['id_1', 'mass_1', 'sink_1_pos_0', 'sink_1_pos_1', 'sink_1_pos_2', \
                            'sink_1_vel_0', 'sink_1_vel_1', 'sink_1_vel_2', \
                            'dm_1', \
                            'dp_1_0', 'dp_1_1', 'dp_1_2', \
                            'dspin_1_0', 'dspin_1_1', 'dspin_1_2', \
                            'f_inst_grav_1_0', 'f_inst_grav_1_1', 'f_inst_grav_1_2', \
                            'f_inst_grav_a_1_0', 'f_inst_grav_a_1_1', 'f_inst_grav_a_1_2', \
                            'f_inst_grav_acc_1_0', 'f_inst_grav_acc_1_1', 'f_inst_grav_acc_1_2', \
                            'dlb_sum_grav_1_0', 'dlb_sum_grav_1_1', 'dlb_sum_grav_1_2', \
                            'dlb_sum_grav_a_1_0', 'dlb_sum_grav_a_1_1', 'dlb_sum_grav_a_1_2', \
                            'dlb_sum_acc_1_0', 'dlb_sum_acc_1_1', 'dlb_sum_acc_1_2', \
                            'depsb_sum_grav_no_mdot_1', \
                            'depsb_sum_grav_acc_no_mdot_1', \
                            'depsb_sum_grav_a_no_mdot_1', \
                            'depsb_sum_acc_no_mdot_1']
            columns_bh_2 = ['id_2', 'mass_2', 'sink_2_pos_0', 'sink_2_pos_1', 'sink_2_pos_2', \
                            'sink_2_vel_0', 'sink_2_vel_1', 'sink_2_vel_2', \
                            'dm_2', \
                            'dp_2_0', 'dp_2_1', 'dp_2_2', \
                            'dp_2_0', 'dp_2_1', 'dp_2_2', \
                            'f_inst_grav_2_0', 'f_inst_grav_2_1', 'f_inst_grav_2_2', \
                            'f_inst_grav_a_2_0', 'f_inst_grav_a_2_1', 'f_inst_grav_a_2_2', \
                            'f_inst_grav_acc_2_0', 'f_inst_grav_acc_2_1', 'f_inst_grav_acc_2_2', \
                            'dlb_grav_sum_2_0', 'dlb_grav_sum_2_1', 'dlb_sum_grav_2_2', \
                            'dlb_sum_grav_a_2_0', 'dlb_sum_grav_a_2_1', 'dlb_sum_grav_a_2_2', \
                            'dlb_sum_acc_2_0', 'dlb_sum_acc_2_1', 'dlb_sum_acc_2_2', \
                            'depsb_sum_grav_no_mdot_2', \
                            'depsb_sum_grav_acc_no_mdot_2', \
                            'depsb_sum_grav_a_no_mdot_2', \
                            'depsb_sum_acc_no_mdot_2']

        torque_array = pd.read_csv(fp + '//%s.txt' %fname, sep=" ", header=None)
        torque_array.columns = columns

        if raw:
            if 'tmin' in kwargs:
                tmin=kwargs['tmin']
            else:
                tmin=min(torque_array['t'])
            if 'tmax' in kwargs:
                tmax = kwargs['tmax']
            else:
                tmax=max(torque_array['t'])

            ind_tmin = torque_array['t']>=tmin
            ind_tmax = torque_array['t']<=tmax
            ind_t = ind_tmin & ind_tmax

            torque_array_return = {}
            for key in columns:
                torque_array_return[key] = torque_array[key][ind_t]
            return(torque_array_return)

        torque_array_ordered = {}
        for torque_array_column in torque_array.columns:
            torque_array_ordered[torque_array_column] = np.zeros(np.shape(torque_array[torque_array_column]))
            torque_array_ordered['t'] = torque_array['t']

        snap_0 = h5.File(fp +  '/snap_000.hdf5')

        bh_ids = []
        if max(snap_0['PartType5']['Masses']) == snap_0['PartType5']['Masses'][0]:
            bh_1_id = snap_0['PartType5']['ParticleIDs'][0]
            bh_1_mass = snap_0['PartType5']['Masses'][0]
            bh_ids.append(bh_1_id)
            bh_2_id = 0
            bh_2_mass = 0
            n_sinks = len(snap_0['PartType5']['ParticleIDs'])
            if n_sinks == 2:
                bh_2_id = snap_0['PartType5']['ParticleIDs'][1]
                bh_2_mass = snap_0['PartType5']['Masses'][1]
        else:
            bh_1_id = snap_0['PartType5']['ParticleIDs'][1]
            bh_1_mass = snap_0['PartType5']['Masses'][1]
            bh_2_id = snap_0['PartType5']['ParticleIDs'][0]
            bh_2_mass = snap_0['PartType5']['Masses'][0]

        for n_sink,id_sink in zip([[1,2],[2,1]],[bh_1_id, bh_2_id]):
            #first turn all sink ids into int, in case they were saved as floats/exp
            torque_array['id_%d' %n_sink[0]] = [int(x) for x in torque_array['id_%d' %n_sink[0]]]
            torque_array['id_%d' %n_sink[1]] = [int(x) for x in torque_array['id_%d' %n_sink[1]]]

            ind_bh_1 = torque_array['id_%d' %n_sink[0]] == id_sink
            ind_bh_2 = torque_array['id_%d' %n_sink[1]] == id_sink

            for column_bh_1, column_bh_2 in zip(columns_bh_1, columns_bh_2):
                torque_array_ordered[column_bh_1][ind_bh_1] = torque_array[column_bh_1][ind_bh_1]
                torque_array_ordered[column_bh_1][ind_bh_2] = torque_array[column_bh_2][ind_bh_2]

                torque_array_ordered[column_bh_2][ind_bh_1] = torque_array[column_bh_2][ind_bh_1]
                torque_array_ordered[column_bh_2][ind_bh_2] = torque_array[column_bh_1][ind_bh_2]

        """ Check if there are discontinuities: when restarting, the time difference can be 0 or -ve """
        inds_to_remove = []
        i_t = 0
        while i_t < (len(torque_array_ordered['t'])-1):
            t_i = torque_array_ordered['t'][i_t]
            t_i_next = torque_array_ordered['t'][i_t+1]
            if (t_i_next-t_i) <= 0:
                index_last_time = i_t
                """ find up until where time is less than or equal to the current index """
                while torque_array_ordered['t'][i_t] <= torque_array_ordered['t'][index_last_time]:
                    inds_to_remove.append(i_t)
                    i_t += 1
                    if i_t >= len(torque_array_ordered['t'])-1:
                        break
            i_t += 1

        """ Now make boolean array for indexing """
        bool_arr = np.ones(np.shape(torque_array_ordered['t']), dtype=bool)
        bool_arr[inds_to_remove] = False
        """ And remove the indeces where duplicates happen """
        for key in torque_array_ordered.keys():
            torque_array_ordered[key] = torque_array_ordered[key][bool_arr]

        """ ACCRETION RATES """
        torque_array_ordered['dm'] = torque_array_ordered['dm_1'] + torque_array_ordered['dm_2']
        torque_array_ordered['mdot_1'] = np.zeros((len(torque_array_ordered['t'])))
        torque_array_ordered['mdot_1'][1:] = torque_array_ordered['dm_1'][1:]/np.diff(torque_array_ordered['t'])
        torque_array_ordered['mdot_2'] = np.zeros((len(torque_array_ordered['t'])))
        torque_array_ordered['mdot_2'][1:] = torque_array_ordered['dm_2'][1:]/np.diff(torque_array_ordered['t'])
        torque_array_ordered['mdot'] = torque_array_ordered['mdot_1'] + torque_array_ordered['mdot_2']

        """ VECTORS CONNECTING SINKS: RELATIVE POSITION AND VELOCITY """
        r_1 = np.array([torque_array_ordered['sink_1_pos_0'], torque_array_ordered['sink_1_pos_1'], torque_array_ordered['sink_1_pos_2']])
        r_2 = np.array([torque_array_ordered['sink_2_pos_0'], torque_array_ordered['sink_2_pos_1'], torque_array_ordered['sink_2_pos_2']])
        torque_array_ordered['r_b'] = r_1 - r_2
        torque_array_ordered['r_b_mag'] = np.linalg.norm(torque_array_ordered['r_b'], axis=0)
        v_1 = np.array([torque_array_ordered['sink_1_vel_0'], torque_array_ordered['sink_1_vel_1'], torque_array_ordered['sink_1_vel_2']])
        v_2 = np.array([torque_array_ordered['sink_2_vel_0'], torque_array_ordered['sink_2_vel_1'], torque_array_ordered['sink_2_vel_2']])
        torque_array_ordered['v_b'] = v_1 - v_2

        """ CORRECT FOR MISSING QUANTITIES IN e=0.200, e=0.400 SIMULATIONS """
        if 'depsb_sum_grav_a_no_mdot_1' not in columns:
            torque_array_ordered['depsb_sum_grav_a_no_mdot_1'] = torque_array_ordered['depsb_sum_grav_a_1'] + torque_array_ordered['dm_1']/torque_array_ordered['r_b_mag']
        if 'depsb_sum_grav_a_no_mdot_2' not in columns:
            torque_array_ordered['depsb_sum_grav_a_no_mdot_2'] = torque_array_ordered['depsb_sum_grav_a_2'] + torque_array_ordered['dm_2']/torque_array_ordered['r_b_mag']
        """ CORRECT FOR MISSING QUANTITIES IN e=0.200, e=0.400 SIMULATIONS """

        if 'tmin' in kwargs:
            tmin=kwargs['tmin']
        else:
            tmin=min(torque_array_ordered['t'])
        if 'tmax' in kwargs:
            tmax = kwargs['tmax']
        else:
            tmax=max(torque_array_ordered['t'])

        ind_tmin = torque_array_ordered['t']>=tmin
        ind_tmax = torque_array_ordered['t']<=tmax
        ind_t = ind_tmin & ind_tmax

        t_orb = misc.T_orbit(a = 1)
        torque_array_ordered['t_orb'] = torque_array_ordered['t']/t_orb

        for key in torque_array_ordered.keys():
            if key == 'r_b' or key == 'v_b':
                for j in [0,1,2]:
                    torque_array_ordered[key][j] = torque_array_ordered[key][j][ind_t]
            else:
                torque_array_ordered[key] = torque_array_ordered[key][ind_t]

        ind_nan = [True]*len(torque_array_ordered['t'])
        for key in torque_array_ordered.keys():
            if key == 'r_b' or key == 'v_b':
                for j in [0,1,2]:
                    ind_nan2 = ~np.isnan(torque_array_ordered[key][j])
                    ind_nan = ind_nan & ind_nan2
            else:
                ind_nan2 = ~np.isnan(torque_array_ordered[key])
                ind_nan = ind_nan & ind_nan2

        for key in torque_array_ordered.keys():
            if key == 'r_b' or key == 'v_b':
                for j in [0,1,2]:
                    torque_array_ordered[key][j] = torque_array_ordered[key][j][ind_nan]
            else:
                torque_array_ordered[key] = torque_array_ordered[key][ind_nan]

        return(torque_array_ordered)


                    


                

            