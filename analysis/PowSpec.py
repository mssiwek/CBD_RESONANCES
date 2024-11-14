""" 
    THIS CLASS TAKES A TIME SERIES, E.G.,
    BINARY EVOLUTION RATES \dot{a}_{\rm b} and \dot{e}_{\rm b},
    AND DOES A FOURIER DECOMPOSITION. 
    IT CAN REMOVE PARTS OF THE POWER SPECTRUM AND 
    CHECK HOW THIS CHANGES THE EVOLUTION RATES.

    created by Magdalena Siwek, last modified 10/2024
"""
import resonances as res
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

class PowSpec:
    def __init__(self, x, y, **kwargs):
        #if time series, x should be in units of binary periods, i.e., t/Pb
        self.x = x
        self.y = y
    
    """ GET THE POWER SPECTRUM """
    def ls(self):
        frequency, power = LombScargle(self.x, self.y, \
                                       normalization='psd').autopower()
        #with normalization='psd', power has units y.unit**2 (e.g., [ebdot]**2 = 1/t**2)
        self.x_ls = frequency
        self.y_ls = power
        return()

    def ls_peaks(self, np=5):
        peaks, _ = find_peaks(self.y_ls, height=0.5)
        #peaks are the indices of the peaks detected in self.y_ls
        #now we want to get the largest peaks (up to rank np):
        y_peaks = self.y_ls[peaks]
        inds_desc = y_peaks.argsort()[::-1] #[::-1] reverses the order
        np_inds = inds_desc[0:np]

        #get x and y vals of np peak indices
        self.x_pks = self.x_ls[np_inds]
        self.y_pks = self.y_ls[np_inds]
        return()


    def nearest_res(self):
        """ FIND THE NEAREST RESONANCES TO A GIVEN RANGE OF FREQUENCIES, 
            SELECT THE ONE WITH THE LOWEST |l-m|.
            IF THERE ARE NO RESONANCES WITHIN \delta f > f_res, 
            RETURN A NaN/SOME OTHER NOTICE FOR THIS RESONANCE """

        


        
