"""
    THIS CLASS TAKES A TIME SERIES, E.G.,
    BINARY EVOLUTION RATES r'\dot{a}_{\rm b} and \dot{e}_{\rm b}',
    AND DOES A FOURIER DECOMPOSITION. 
    IT CAN REMOVE PARTS OF THE POWER SPECTRUM AND 
    CHECK HOW THIS CHANGES THE EVOLUTION RATES.

    created by Magdalena Siwek, last modified 10/2024
"""
import resonances as res
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import numpy as np
import misc
import plot_utils as put

class PowSpec:
    def __init__(self, x, y, **kwargs):
        #if time series, x should be in units of binary periods, i.e., t/Pb
        self.x = x
        self.y = y
    
    """ GET THE POWER SPECTRUM """
    def ls(self, xlim=None):
        frequency, power = LombScargle(self.x, self.y, \
                                       normalization='psd').autopower()
        
        # power = np.sqrt(power * len(self.y)) #convert into fft units
        #with normalization='psd', power has units y.unit**2 (e.g., [ebdot]**2 = 1/t**2)
        
        self.x_ls = frequency
        self.y_ls = power
        if xlim is None:
            pass
        else:
            inds_lower = self.x_ls >= xlim[0]
            inds_upper = self.x_ls <= xlim[1]
            inds = inds_upper & inds_lower 
            self.x_ls = self.x_ls[inds]
            self.y_ls = self.y_ls[inds]

        return()

    def ls_peaks(self, npeaks=5):
        peaks, _ = find_peaks(self.y_ls)
        #peaks are the indices of the peaks detected in self.y_ls
        #now we want to get the largest peaks (up to rank npeaks):
        y_peaks = self.y_ls[peaks]
        x_peaks = self.x_ls[peaks]
        inds_desc = y_peaks.argsort()[::-1] #[::-1] reverses the order
        np_inds = inds_desc[0:npeaks]

        #get x and y vals of npeaks indices
        self.x_pks = x_peaks[np_inds]
        self.y_pks = y_peaks[np_inds]
        return()


    def nearest_res(self, lim_m=10, lim_l=10, width=5, \
                    res_types=['CR', 'ILR', 'OLR']):
        """ FIND THE NEAREST RESONANCES TO A GIVEN RANGE OF FREQUENCIES,
            WITHIN THE RANGE 2*WIDTH*FREQUENCY AROUND EACH FREQUENCY. 
            SELECT THE ONE WITH THE LOWEST |l-m|.
            IF THERE ARE NO RESONANCES WITHIN delta f > f_res, 
            RETURN A NaN/SOME OTHER NOTICE FOR THIS RESONANCE """

        #GET ALL THE FREQUENCIES ASSOCIATED WITH CR, ILR AND OLR
        res_dict = {}
        res_dict['res_omega'] = []
        res_dict['lm'] = []
        res_dict['res_type'] = []
        for res_type in res_types:
            for l in range(0,lim_l):
                for m in range(0,lim_m):
                    r = eval('res.r_%s(l,m)' %res_type)
                    if r == np.nan:
                        continue
                    else:
                        res_dict['res_omega'].append(res.omega(r))
                        res_dict['lm'].append((l,m))
                        res_dict['res_type'].append(res_type)
        
        #NOW LOOP OVER PEAK FREQUENCIES AND FIND THE NEAREST RESONANCES
        self.match_res = {}
        keywords = ['n', 'omega', 'res_type', 'res_omega', 'l', 'm']
        for keyword in keywords:
            self.match_res[keyword] = []

        for n,omega in enumerate(self.x_pks):
            
            ind_lower = res_dict['res_omega'] >= (omega - omega*width)
            ind_upper = res_dict['res_omega'] <= (omega + omega*width)
            ind_omega_range = ind_lower & ind_upper

            # There are resonances within the specified window
            if np.any(ind_omega_range):
                #find the one with the smallest |l-m|
                all_lm = np.array(res_dict['lm'])[ind_omega_range]
                diff_lm = all_lm[:,0] - all_lm[:,1]
                #take the absolute value
                diff_lm[diff_lm<0] *= -1
                #find the resonance with the minimum value |l-m|
                ind_min_lm = diff_lm == min(diff_lm)

                # Pick the most relevant resonance to return for current peak
                if ind_min_lm.sum() > 1:
                    # finding more than one resonance with minimum |l-m|, 
                    # pick the one with frequency closest to omega
                    res_omegas = np.array(res_dict['res_omega'])[ind_omega_range][ind_min_lm]
                    ind_final = misc.closest_idx(res_omegas, omega)
                else:
                    # only one resonance in the vicinity is found
                    ind_final = np.where(ind_min_lm)[0][0]
                
                self.match_res['n'].append(n)
                self.match_res['omega'].append(omega)
                self.match_res['res_type'].append(np.array(res_dict['res_type'])[ind_omega_range][ind_final])
                self.match_res['res_omega'].append(np.array(res_dict['res_omega'])[ind_omega_range][ind_final])
                self.match_res['l'].append(np.array(res_dict['lm'])[ind_omega_range][:,0][ind_final])
                self.match_res['m'].append(np.array(res_dict['lm'])[ind_omega_range][:,1][ind_final])

            # No resonances within the specified window
            else:
                continue
        
        for key in self.match_res.keys():
            self.match_res[key] = np.array(self.match_res[key])
        return()

    # NOW ASSIGN PLOTTING COLORS FOR EACH OF THE RESONANCES
    def get_colors(self, base_colors=None):
        if not hasattr(self, 'match_res'):
            exit("PowSpec must have attribute 'match_res' to assign colors.")
        
        if base_colors is not None:
            self.base_colors = base_colors
        else:
            self.base_colors = {
                'ILR': 'seagreen',\
                'OLR': 'royalblue',\
                'CR': 'mediumvioletred',\
            }
            
        
        self.match_res['colors'] = np.empty([len(self.match_res['res_type']),3])
        for key in self.base_colors.keys():
            ind_res = np.where(self.match_res['res_type'] == key)[0]
            if len(ind_res) > 0:
                #we have found this resonance
                cfacs = np.linspace(1.2,0.3,len(ind_res))
                for cfac,ind in zip(cfacs,ind_res):
                    self.match_res['colors'][ind] = \
                        put.lighten_color(self.base_colors[key], cfac)
            else:
                continue
        return()