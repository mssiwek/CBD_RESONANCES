""" 
    THIS CLASS TAKES A TIME SERIES, E.G.,
    BINARY EVOLUTION RATES \dot{a}_{\rm b} and \dot{e}_{\rm b},
    AND DOES A FOURIER DECOMPOSITION. 
    IT CAN THEN FIT THE TIME SERIES WITH FOURIER COMPONENTS,
    AND REMOVE COMPONENTS TO CHECK IF THESE CONTRIBUTE SIGNIFICANTLY
    TO THE SIGN OF THE EVOLUTION RATES.

    created by Magdalena Siwek, last modified 10/2024
"""

class PowDecomp:
    def __init__(self, x, y, **kwargs):
        self.x = x
        self.y = y
    
    """ GET THE POWER SPECTRUM """
    def powspec(self):

    
    """ FIT THE POWER SPECTRUM """
    def fit(x, y):
        