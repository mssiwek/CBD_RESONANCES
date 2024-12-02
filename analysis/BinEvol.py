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
import torques_eb_lb as tq
import misc
import sys
import time
import h5py as h5

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
                    - self.SimInit.depsb/self.SimInit.epsb \
                    - (2. * self.SimInit.dlb/self.SimInit.lb))
        else:
            ebdot = (1 - self.eb**2) \
                    * (2. * self.SimInit.mbdot/self.mb \
                    - self.SimInit.depsb/self.SimInit.epsb \
                    - (2. * self.SimInit.dlb/self.SimInit.lb))
        return(ebdot)

    def abdot(self):
        abdot = self.SimInit.mbdot/self.mb \
                - self.SimInit.depsb/self.SimInit.epsb
        return(abdot) 


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
        
    
    def evol_from_txt(self):
        """ THIS MEANS WE ARE SUPPLYING A FILEPATH VIA ext,
            AND GET dlb AND depsb FROM torques_eb_lb FILE """
        
        """ LOAD SIMULATION OUTPUT AND PARAMETERS """
        """ ------------------------------------- """
        """ ------------------------------------- """
        start_time = time.time()
        tqfl = tq.read_torques_eb_lb_txt(self.ext, orbit_averaged=True, n_orbit_average = 1, \
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



        """ --------------------------------------- """
        """ --------------------------------------- """
        """ GETTING ANGULAR MOMENTUM RATE OF CHANGE """
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
        """ GETTING ANGULAR MOMENTUM RATE OF CHANGE """
        """ --------------------------------------- """
        """ --------------------------------------- """



        """ ----------------------------- """
        """ ----------------------------- """
        """ GETTING ENERGY RATE OF CHANGE """
        if not self.fg_cav:
            """ f_ext = f_grav ONLY """
            self.depsb1 = tqfl['depsb_sum_grav_no_mdot_1']
            self.depsb2 = tqfl['depsb_sum_grav_no_mdot_2']
        else:
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
        """ GETTING ENERGY RATE OF CHANGE """
        """ ----------------------------- """
        """ ----------------------------- """


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
        
        """ CALCULATE TOTAL RATE OF CHANGE OF 
            ANGULAR MOMENTUM AND ENERGY """
        """ --------------------------------- """
        """ --------------------------------- """
        self.dlb = s1 * self.dlb1 + s2 * self.dlb2
        self.depsb = s1 * self.depsb1 + s2 * self.depsb2
        """ --------------------------------- """
        """ --------------------------------- """
        """ CALCULATE TOTAL RATE OF CHANGE OF 
            ANGULAR MOMENTUM AND ENERGY """
        
        if not self.fg_cav:
            if self.acc:
                """ NOW WE ARE ALSO ADDING EFFECT OF ACCRETION TO THE 
                    RATE OF CHANGE OF EPSILON (BINARY SPECIFIC ENERGY) """
                self.depsb += - tqfl['dm']/tqfl['r_b_mag']
                self.mbdot = tqfl['dm']
            else:
                self.mbdot = 0
        else:
            self.mbdot = 0
        print("remainder of evol_from_txt took %s seconds ---" % (time.time() - start_time))

        return(self)
        

    def evol_from_snap(self, sn):
        print("Note: ebdot and abdot maps do not include accretion contributions.")
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
        self.dlb = rb[0]*self.fgrav[1] - rb[1]*self.fgrav[0]
        """ GETTING ANGULAR MOMENTUM RATE OF CHANGE """
        """ --------------------------------------- """
        """ --------------------------------------- """
    

        """ ----------------------------- """
        """ ----------------------------- """
        """ GETTING ENERGY RATE OF CHANGE """
        vb = sn['PartType5']['Velocities'][0] - sn['PartType5']['Velocities'][1]
        self.depsb = vb[0]*self.fgrav[0]+vb[1]*self.fgrav[1]
        """ GETTING ENERGY RATE OF CHANGE """
        """ ----------------------------- """
        """ ----------------------------- """

        self.mbdot = 0
        
        return(self)

                


            

        